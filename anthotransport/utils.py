from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import shlex
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Sequence


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def atomic_write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", dir=path.parent, delete=False
    ) as handle:
        handle.write(text)
        temporary = Path(handle.name)
    os.replace(temporary, path)


def atomic_write_json(path: Path, value: Any) -> None:
    atomic_write_text(path, json.dumps(value, indent=2, sort_keys=True) + "\n")


def atomic_write_tsv(
    path: Path, rows: Iterable[dict[str, Any]], fieldnames: Sequence[str]
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w", encoding="utf-8", newline="", dir=path.parent, delete=False
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
        temporary = Path(handle.name)
    os.replace(temporary, path)


def read_tsv(path: Path, required: Iterable[str] = ()) -> list[dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(f"TSV file not found: {path}")
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = set(required) - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                f"{path}: missing columns: {', '.join(sorted(missing))}"
            )
        return list(reader)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def fingerprint_path(path: Path) -> dict[str, Any]:
    """Fingerprint files by content and directories by a deterministic inventory."""
    path = path.resolve()
    if not path.exists():
        return {"path": str(path), "missing": True}
    if path.is_file():
        stat = path.stat()
        return {
            "path": str(path),
            "kind": "file",
            "size": stat.st_size,
            "sha256": sha256_file(path),
        }
    entries = []
    for child in sorted(item for item in path.rglob("*") if item.is_file()):
        stat = child.stat()
        entries.append(
            {
                "path": str(child.relative_to(path)),
                "size": stat.st_size,
                "mtime_ns": stat.st_mtime_ns,
            }
        )
    encoded = json.dumps(entries, sort_keys=True, separators=(",", ":")).encode()
    return {
        "path": str(path),
        "kind": "directory",
        "files": len(entries),
        "inventory_sha256": hashlib.sha256(encoded).hexdigest(),
    }


def signature(value: Any) -> str:
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def fasta_ids(path: Path) -> list[str]:
    ids: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    return ids


def require_unique(values: Sequence[str], label: str) -> None:
    seen: set[str] = set()
    duplicates: set[str] = set()
    for value in values:
        if value in seen:
            duplicates.add(value)
        seen.add(value)
    if duplicates:
        preview = ", ".join(sorted(duplicates)[:20])
        raise ValueError(f"Duplicate IDs in {label}: {preview}")


def resolve_tree_id(leaf_names: set[str], node_id: str) -> str | None:
    candidates = [node_id]
    sanitized = node_id.replace("@", "_")
    if sanitized != node_id:
        candidates.append(sanitized)
    for candidate in candidates:
        if candidate in leaf_names:
            return candidate
    return None


def tree_ids_from_newick_text(tree: str) -> set[str]:
    # Used only for lightweight checks. Tree algorithms use Bio.Phylo.
    return set(re.findall(r"(?:^|[(,])'?([^'():,;]+)'?(?=[:),])", tree))


def format_command(command: Sequence[str]) -> str:
    return shlex.join(str(part) for part in command)


class CommandRunner:
    def __init__(self, log_file: Path, dry_run: bool = False):
        self.log_file = log_file
        self.dry_run = dry_run

    def run(
        self,
        command: Sequence[str],
        *,
        stdout_path: Path | None = None,
        stderr_path: Path | None = None,
        cwd: Path | None = None,
    ) -> None:
        rendered = format_command(command)
        self.log_file.parent.mkdir(parents=True, exist_ok=True)
        with self.log_file.open("a", encoding="utf-8") as log:
            log.write(f"[{utc_now()}] cwd={cwd or Path.cwd()} command={rendered}\n")
        print(f"$ {rendered}")
        if self.dry_run:
            return

        stdout_handle = stdout_path.open("wb") if stdout_path else None
        same_log = stdout_path is not None and stderr_path == stdout_path
        stderr_handle = None if same_log else (stderr_path.open("wb") if stderr_path else None)
        try:
            subprocess.run(
                [str(part) for part in command],
                cwd=cwd,
                stdout=stdout_handle,
                stderr=subprocess.STDOUT if same_log else stderr_handle,
                check=True,
            )
        finally:
            if stdout_handle:
                stdout_handle.close()
            if stderr_handle:
                stderr_handle.close()


def executable_version(command: Sequence[str]) -> str:
    try:
        result = subprocess.run(
            list(command), text=True, capture_output=True, timeout=20, check=False
        )
    except (OSError, subprocess.SubprocessError) as exc:
        return f"unavailable: {exc}"
    combined = (result.stdout + "\n" + result.stderr).strip()
    return combined.splitlines()[0] if combined else f"exit {result.returncode}"

