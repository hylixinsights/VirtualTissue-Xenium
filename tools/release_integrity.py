"""Bind local release checks and a reviewable GitHub payload to exact source bytes.

snapshot(repo) returns a canonical inventory of Git-tracked files (including
working-tree deletions) and non-ignored untracked files. Git executable modes are
included. Symlinks, submodules and special files are refused, not dereferenced.
Ignored environments, outputs and secrets are outside this source inventory.

compare_snapshots(before, after) raises IntegrityError on any changed inventory,
content, executable mode or HEAD. write_github_payload(repo, output_dir,
tested_snapshot, base_sha) requires unchanged source and HEAD == base_sha. It
creates candidate-source.tar.gz, candidate.patch, source-snapshot.json,
candidate-metadata.json and SHA256SUMS in a new directory. The patch includes
untracked files and deletions and is applied to a fresh temporary Git index to
verify its exact resulting tree. Neither the real index nor source files nor
repository object database are changed. There is no commit, push or publication.

This is source integrity, not a claim that tests ran or that biology is valid.
The caller must revoke any previous green status before running its checks,
record test/environment evidence, and bind that result to tree_sha256 and the
returned artifact hashes. Executable mode means Git's owner-executable flag;
ownership, ACLs and other local permission bits are not release source content.
"""
from __future__ import annotations

from contextlib import contextmanager
import gzip
import hashlib
import io
import json
import os
from pathlib import Path, PurePosixPath
import re
import stat
import subprocess
import tarfile
import tempfile
import tomllib

VERSION = "release-integrity-1"


class IntegrityError(RuntimeError):
    """A release no longer matches its declared source or base."""


def _canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True, allow_nan=False).encode()


def _sha256(data):
    return hashlib.sha256(data).hexdigest()


def _git(repo, *arguments, env=None, input=None):
    process_env = os.environ.copy()
    # Do not accidentally address another worktree or a caller's temporary index.
    for key in ("GIT_DIR", "GIT_WORK_TREE", "GIT_INDEX_FILE", "GIT_COMMON_DIR",
                "GIT_OBJECT_DIRECTORY", "GIT_ALTERNATE_OBJECT_DIRECTORIES"):
        process_env.pop(key, None)
    process_env.update(GIT_OPTIONAL_LOCKS="0", LC_ALL="C")
    process_env.update(env or {})
    command = ["git", "-C", str(repo), "-c", "core.fsmonitor=false", "-c", "core.filemode=true",
               "-c", "core.autocrlf=false", "-c", "core.safecrlf=false", *arguments]
    result = subprocess.run(command, input=input, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=process_env)
    if result.returncode:
        raise IntegrityError(f"Git integrity check failed ({arguments[0]}): {result.stderr.decode(errors='replace').strip()}")
    return result.stdout


def _root(repo):
    root = Path(repo).resolve(strict=True)
    actual = Path(os.fsdecode(_git(root, "rev-parse", "--show-toplevel")).strip()).resolve()
    if actual != root:
        raise IntegrityError("Pass the repository root, not a parent or subdirectory.")
    return root


def _safe_path(value):
    if not isinstance(value, str) or not value or "\x00" in value or "\\" in value:
        raise IntegrityError("Source paths must be safe relative POSIX paths.")
    path = PurePosixPath(value)
    if (path.is_absolute() or path.as_posix() != value or any(part in ("", ".", "..") for part in path.parts)
            or any(part.casefold() == ".git" for part in path.parts) or re.match(r"^[A-Za-z]:", value)):
        raise IntegrityError(f"Unsafe source path: {value!r}.")
    try:
        value.encode("utf-8")
    except UnicodeEncodeError as exc:
        raise IntegrityError("Source filenames must be valid UTF-8.") from exc
    return path


def _check_parents(root, relative):
    path = root
    for part in _safe_path(relative).parts:
        path = path / part
        if path.is_symlink():
            raise IntegrityError(f"Symlinks are not allowed in release source: {relative!r}.")
    return path


def _stat_identity(value):
    return (value.st_dev, value.st_ino, value.st_size, value.st_mtime_ns, value.st_ctime_ns, value.st_mode)


def _file_record(root, relative):
    path = _check_parents(root, relative)
    try:
        initial = path.lstat()
    except FileNotFoundError:
        return {"status": "deleted", "sha256": None, "mode": None, "size": None}
    if not stat.S_ISREG(initial.st_mode):
        raise IntegrityError(f"Only regular source files are supported: {relative!r}.")
    try:
        fd = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0))
    except FileNotFoundError:
        return {"status": "deleted", "sha256": None, "mode": None, "size": None}
    with os.fdopen(fd, "rb") as handle:
        before = os.fstat(handle.fileno())
        if not stat.S_ISREG(before.st_mode):
            raise IntegrityError(f"Only regular source files are supported: {relative!r}.")
        checksum = hashlib.sha256()
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            checksum.update(chunk)
        after = os.fstat(handle.fileno())
    _check_parents(root, relative)
    try:
        current = path.stat()
    except FileNotFoundError as exc:
        raise IntegrityError(f"Source changed while hashing {relative!r}.") from exc
    if _stat_identity(before) != _stat_identity(after) or _stat_identity(before) != _stat_identity(current):
        raise IntegrityError(f"Source changed while hashing {relative!r}.")
    return {"status": "file", "sha256": checksum.hexdigest(), "size": before.st_size,
            "mode": "100755" if before.st_mode & stat.S_IXUSR else "100644"}


def snapshot(repo):
    """Return {version, head, files, tree_sha256}; no Git/index/source mutation.

    files maps relative paths to {status, sha256, size, mode}. Missing tracked
    files have status='deleted' and null content/mode. tree_sha256 hashes the
    canonical files mapping; head is compared independently by compare_snapshots.
    """
    root = _root(repo)
    head = _git(root, "rev-parse", "HEAD").decode().strip()
    listed = _git(root, "ls-files", "--cached", "--others", "--exclude-standard", "-z")
    names = sorted({os.fsdecode(part) for part in listed.split(b"\0") if part})
    files = {name: _file_record(root, name) for name in names}
    listed_after = _git(root, "ls-files", "--cached", "--others", "--exclude-standard", "-z")
    if set(listed.split(b"\0")) != set(listed_after.split(b"\0")):
        raise IntegrityError("Source paths changed while source was being inventoried.")
    if _git(root, "rev-parse", "HEAD").decode().strip() != head:
        raise IntegrityError("HEAD changed while source was being inventoried.")
    return {"version": VERSION, "head": head, "files": files, "tree_sha256": _sha256(_canonical(files))}


def _validate_snapshot(value):
    if not isinstance(value, dict) or value.get("version") != VERSION or not isinstance(value.get("files"), dict):
        raise IntegrityError("Expected a release-integrity source snapshot.")
    if value.get("tree_sha256") != _sha256(_canonical(value["files"])):
        raise IntegrityError("Snapshot tree hash does not match its recorded file inventory.")
    if not isinstance(value.get("head"), str) or not re.fullmatch(r"[0-9a-f]{40}|[0-9a-f]{64}", value["head"]):
        raise IntegrityError("Snapshot HEAD must be a full Git commit ID.")
    for name, record in value["files"].items():
        _safe_path(name)
        if not isinstance(record, dict) or set(record) != {"status", "sha256", "mode", "size"}:
            raise IntegrityError(f"Invalid source record for {name!r}.")
        if record["status"] == "deleted":
            if any(record[key] is not None for key in ("sha256", "mode", "size")):
                raise IntegrityError("Deleted source records cannot contain current file data.")
        elif (record["status"] != "file" or record["mode"] not in ("100644", "100755")
              or type(record["size"]) is not int or record["size"] < 0
              or not isinstance(record["sha256"], str) or not re.fullmatch(r"[0-9a-f]{64}", record["sha256"])):
            raise IntegrityError(f"Invalid current source record for {name!r}.")


def compare_snapshots(old, new):
    """Raise IntegrityError on changed source or HEAD; return True on identity."""
    _validate_snapshot(old)
    _validate_snapshot(new)
    if old["head"] != new["head"]:
        raise IntegrityError("Repository HEAD changed after the tested source was captured.")
    changes = [path for path in sorted(set(old["files"]) | set(new["files"]))
               if old["files"].get(path) != new["files"].get(path)]
    if changes:
        preview = ", ".join(repr(path) for path in changes[:12])
        raise IntegrityError(f"Source changed after snapshot ({len(changes)} paths): {preview}.")
    return True


@contextmanager
def _isolated_index(root):
    """Temporary index AND object store; read existing objects via alternates."""
    raw_objects = Path(os.fsdecode(_git(root, "rev-parse", "--git-path", "objects")).strip())
    existing = (root / raw_objects).resolve() if not raw_objects.is_absolute() else raw_objects.resolve()
    with tempfile.TemporaryDirectory(prefix="vtx-release-index-") as directory:
        temporary = Path(directory)
        objects = temporary / "objects"
        objects.mkdir()
        yield temporary, {"GIT_INDEX_FILE": str(temporary / "candidate.index"),
            "GIT_OBJECT_DIRECTORY": str(objects),
            "GIT_ALTERNATE_OBJECT_DIRECTORIES": json.dumps(str(existing))}


def _index_records(root, env):
    result = {}
    for entry in _git(root, "ls-files", "--stage", "-z", env=env).split(b"\0"):
        if not entry:
            continue
        fields, name = entry.split(b"\t", 1)
        mode, oid, stage = fields.decode("ascii").split()
        relative = os.fsdecode(name)
        _safe_path(relative)
        if stage != "0" or mode not in ("100644", "100755") or relative in result:
            raise IntegrityError("The candidate index contains unmerged or unsupported objects.")
        result[relative] = {"mode": mode, "oid": oid}
    return result


def _write_json(path, value):
    with path.open("x", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True, ensure_ascii=True, allow_nan=False)
        handle.write("\n")


def _artifact(path):
    checksum = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            checksum.update(chunk)
    return {"path": path.name, "sha256": checksum.hexdigest(), "size": path.stat().st_size}


def write_github_payload(repo, output_dir, tested_snapshot, base_sha):
    """Create a new source tar, exact-base binary patch and checksum metadata.

    Return metadata plus metadata_file and checksums_file artifact records. All
    source artifacts are tied to tested_snapshot.tree_sha256. A partial output
    directory after any error is not valid: candidate-metadata.json is written
    only after all checks pass, and callers must require successful return.
    """
    root = _root(repo)
    _validate_snapshot(tested_snapshot)
    if not isinstance(base_sha, str) or not re.fullmatch(r"[0-9a-f]{40}|[0-9a-f]{64}", base_sha):
        raise IntegrityError("base_sha must be an exact full commit ID.")
    compare_snapshots(tested_snapshot, snapshot(root))
    if tested_snapshot["head"] != base_sha:
        raise IntegrityError("The requested base commit is not the current repository HEAD.")
    output = Path(output_dir).absolute()
    # A payload inside the candidate would itself change the tested inventory.
    if output.resolve().is_relative_to(root):
        raise IntegrityError("Write the GitHub payload outside the source repository.")
    output.mkdir(parents=True, exist_ok=False)
    expected = {name: row for name, row in tested_snapshot["files"].items() if row["status"] == "file"}
    with _isolated_index(root) as (temporary, env):
        _git(root, "read-tree", base_sha, env=env)
        _git(root, "add", "-A", "--", ".", env=env)
        indexed = _index_records(root, env)
        if set(indexed) != set(expected):
            raise IntegrityError("Candidate index and tested source paths differ; check staged ignored files or changed ignore rules.")
        candidate_git_tree = _git(root, "write-tree", env=env).decode().strip()
        patch = _git(root, "diff", "--cached", "--binary", "--full-index", "--no-ext-diff", "--no-textconv",
                     "--no-renames", "--src-prefix=a/", "--dst-prefix=b/", base_sha, "--", env=env)
        (output / "candidate.patch").write_bytes(patch)
        verification_env = dict(env, GIT_INDEX_FILE=str(temporary / "verification.index"))
        _git(root, "read-tree", base_sha, env=verification_env)
        if patch:
            _git(root, "apply", "--check", "--cached", "--binary", str(output / "candidate.patch"), env=verification_env)
            _git(root, "apply", "--cached", "--binary", str(output / "candidate.patch"), env=verification_env)
        applied_tree = _git(root, "write-tree", env=verification_env).decode().strip()
        if applied_tree != candidate_git_tree:
            raise IntegrityError("Patch verification did not reproduce the candidate Git tree.")

        version = None
        with (output / "candidate-source.tar.gz").open("xb") as raw_archive:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw_archive, mtime=0) as compressed:
                with tarfile.open(fileobj=compressed, mode="w", format=tarfile.PAX_FORMAT) as archive:
                    for relative, record in sorted(indexed.items()):
                        data = _git(root, "cat-file", "blob", record["oid"], env=env)
                        wanted = expected[relative]
                        if (record["mode"] != wanted["mode"] or len(data) != wanted["size"]
                                or _sha256(data) != wanted["sha256"]):
                            raise IntegrityError(f"Staged bytes/mode differ from tested source: {relative!r}; check Git filters or source changes.")
                        if relative == "pyproject.toml":
                            version = tomllib.loads(data.decode("utf-8")).get("project", {}).get("version")
                        info = tarfile.TarInfo("candidate/" + relative)
                        info.size, info.mode, info.mtime = len(data), int(record["mode"][-3:], 8), 0
                        info.uid = info.gid = 0
                        info.uname = info.gname = ""
                        archive.addfile(info, io.BytesIO(data))
    compare_snapshots(tested_snapshot, snapshot(root))
    _write_json(output / "source-snapshot.json", tested_snapshot)
    artifacts = {name: _artifact(output / name) for name in
                 ("candidate-source.tar.gz", "candidate.patch", "source-snapshot.json")}
    metadata = dict(version=VERSION, base_sha=base_sha, source_tree_sha256=tested_snapshot["tree_sha256"],
        candidate_git_tree=candidate_git_tree, project_version=version, archive_root="candidate",
        patch_verification="applied_to_exact_base_in_isolated_index_and_matched_candidate_tree",
        source_scope="Git tracked plus non-ignored untracked files; regular files and Git executable mode; tracked deletions recorded",
        publication_status="local_review_payload_only_not_committed_or_pushed", artifacts=artifacts)
    _write_json(output / "candidate-metadata.json", metadata)
    metadata_file = _artifact(output / "candidate-metadata.json")
    with (output / "SHA256SUMS").open("x", encoding="ascii") as handle:
        for record in [*artifacts.values(), metadata_file]:
            handle.write(f"{record['sha256']}  {record['path']}\n")
    return dict(metadata, metadata_file=metadata_file, checksums_file=_artifact(output / "SHA256SUMS"))
