"""Release payloads must contain exactly the tested source, including new files."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import tarfile

import pytest

SPEC = importlib.util.spec_from_file_location("release_integrity", Path(__file__).resolve().parents[1] / "tools/release_integrity.py")
integrity = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(integrity)


def git(repo, *args, env=None):
    return subprocess.check_output(["git", "-C", str(repo), *args], env=env, stderr=subprocess.PIPE)


@pytest.fixture
def repo(tmp_path):
    value = tmp_path / "candidate"
    value.mkdir()
    git(value, "init", "-q")
    git(value, "config", "user.name", "Integrity fixture")
    git(value, "config", "user.email", "fixture@example.invalid")
    git(value, "config", "commit.gpgsign", "false")
    (value / "tracked.txt").write_text("base content\n")
    (value / "delete.txt").write_text("remove this from candidate\n")
    (value / "binary.dat").write_bytes(b"\x00base\xff\x01")
    (value / "pyproject.toml").write_text('[project]\nname="fixture"\nversion="3.0.0a3"\n')
    (value / ".gitignore").write_text("ignored/\n.env\n")
    git(value, "add", ".")
    git(value, "commit", "-qm", "base fixture")
    return value


@pytest.mark.parametrize("change", ["content", "mode", "addition", "deletion"])
def test_source_inventory_detects_every_release_change(repo, change):
    before = integrity.snapshot(repo)
    target = repo / "tracked.txt"
    if change == "content":
        target.write_text("changed content\n")
    elif change == "mode":
        target.chmod(0o755)
    elif change == "addition":
        (repo / "new file.txt").write_text("untracked source\n")
    else:
        target.unlink()
    after = integrity.snapshot(repo)
    if change == "deletion":
        assert after["files"]["tracked.txt"] == {"status": "deleted", "sha256": None, "mode": None, "size": None}
    if change == "mode":
        assert after["files"]["tracked.txt"]["mode"] == "100755"
    with pytest.raises(integrity.IntegrityError, match="Source changed"):
        integrity.compare_snapshots(before, after)


def test_snapshot_is_stable_and_rejects_tampered_inventory(repo):
    before = integrity.snapshot(repo)
    assert integrity.compare_snapshots(before, integrity.snapshot(repo)) is True
    before["files"]["tracked.txt"]["sha256"] = "0" * 64
    with pytest.raises(integrity.IntegrityError, match="tree hash"):
        integrity.compare_snapshots(before, integrity.snapshot(repo))


def test_symlink_file_and_symlink_parent_are_never_followed(repo, tmp_path):
    outside = tmp_path / "outside.txt"
    outside.write_text("do not package outside source\n")
    (repo / "link.txt").symlink_to(outside)
    with pytest.raises(integrity.IntegrityError, match="Symlinks"):
        integrity.snapshot(repo)
    (repo / "link.txt").unlink()
    nested = repo / "nested"
    nested.mkdir()
    (nested / "file.txt").write_text("tracked nested file\n")
    git(repo, "add", "nested/file.txt")
    (nested / "file.txt").unlink()
    nested.rmdir()
    nested.symlink_to(tmp_path, target_is_directory=True)
    with pytest.raises(integrity.IntegrityError, match="Symlinks"):
        integrity.snapshot(repo)


def test_payload_includes_untracked_binary_mode_and_deletion_without_real_git_changes(repo, tmp_path):
    base = git(repo, "rev-parse", "HEAD").decode().strip()
    (repo / "tracked.txt").write_text("candidate content\n")
    (repo / "binary.dat").write_bytes(b"\x00candidate\xff\x03\x00")
    (repo / "delete.txt").unlink()
    (repo / "new script.sh").write_text("#!/bin/sh\nexit 0\n")
    (repo / "new script.sh").chmod(0o755)
    (repo / ".env").write_text("not included\n")
    source = integrity.snapshot(repo)
    index_before = (repo / ".git/index").read_bytes()
    objects_before = {str(p.relative_to(repo / ".git/objects")): p.read_bytes()
                      for p in (repo / ".git/objects").rglob("*") if p.is_file()}
    output = tmp_path / "handoff"
    metadata = integrity.write_github_payload(repo, output, source, base)
    assert metadata["source_tree_sha256"] == source["tree_sha256"]
    assert metadata["base_sha"] == base and metadata["project_version"] == "3.0.0a3"
    assert metadata["publication_status"] == "local_review_payload_only_not_committed_or_pushed"
    assert (repo / ".git/index").read_bytes() == index_before
    objects_after = {str(p.relative_to(repo / ".git/objects")): p.read_bytes()
                     for p in (repo / ".git/objects").rglob("*") if p.is_file()}
    assert objects_after == objects_before
    assert integrity.compare_snapshots(source, integrity.snapshot(repo))
    patch = (output / "candidate.patch").read_bytes()
    assert b"new script.sh" in patch and b"deleted file mode" in patch and b"GIT binary patch" in patch
    with tarfile.open(output / "candidate-source.tar.gz", "r:gz") as archive:
        names = set(archive.getnames())
        assert "candidate/new script.sh" in names
        assert "candidate/delete.txt" not in names and "candidate/.env" not in names
        assert archive.getmember("candidate/new script.sh").mode == 0o755
        for name, record in source["files"].items():
            if record["status"] == "file":
                assert hashlib.sha256(archive.extractfile("candidate/" + name).read()).hexdigest() == record["sha256"]
    assert json.loads((output / "source-snapshot.json").read_text()) == source
    for line in (output / "SHA256SUMS").read_text().splitlines():
        checksum, name = line.split("  ")
        assert hashlib.sha256((output / name).read_bytes()).hexdigest() == checksum
    # Independently verify applicability to the exact base using another index.
    env = dict(os.environ, GIT_INDEX_FILE=str(tmp_path / "external-check.index"))
    git(repo, "read-tree", base, env=env)
    git(repo, "apply", "--check", "--cached", "--binary", str(output / "candidate.patch"), env=env)


def test_stale_source_wrong_base_and_existing_output_refuse_payload(repo, tmp_path):
    source = integrity.snapshot(repo)
    with pytest.raises(integrity.IntegrityError, match="base commit"):
        integrity.write_github_payload(repo, tmp_path / "wrong-base", source, "0" * 40)
    assert not (tmp_path / "wrong-base").exists()
    with pytest.raises(integrity.IntegrityError, match="outside"):
        integrity.write_github_payload(repo, repo / "payload", source, source["head"])
    directory = tmp_path / "existing"
    directory.mkdir()
    (directory / "retained.txt").write_text("do not overwrite\n")
    with pytest.raises(FileExistsError):
        integrity.write_github_payload(repo, directory, source, source["head"])
    assert (directory / "retained.txt").read_text() == "do not overwrite\n"
    (repo / "tracked.txt").write_text("changed after testing\n")
    with pytest.raises(integrity.IntegrityError, match="Source changed"):
        integrity.write_github_payload(repo, tmp_path / "stale", source, source["head"])
    assert not (tmp_path / "stale").exists()


def test_change_during_packaging_never_gets_success_metadata(repo, tmp_path, monkeypatch):
    source = integrity.snapshot(repo)
    original_git = integrity._git
    changed = False
    def changing_git(root, *arguments, **kwargs):
        nonlocal changed
        if arguments[0] == "cat-file" and not changed:
            (repo / "tracked.txt").write_text("changed while packaging\n")
            changed = True
        return original_git(root, *arguments, **kwargs)
    monkeypatch.setattr(integrity, "_git", changing_git)
    with pytest.raises(integrity.IntegrityError, match="Source changed"):
        integrity.write_github_payload(repo, tmp_path / "changed", source, source["head"])
    assert not (tmp_path / "changed/candidate-metadata.json").exists()


def test_git_clean_filter_changes_cannot_silently_diverge_tar_and_tested_source(repo, tmp_path):
    (repo / ".gitattributes").write_text("filtered.txt text eol=lf\n")
    (repo / "filtered.txt").write_bytes(b"source has CRLF\r\n")
    source = integrity.snapshot(repo)
    with pytest.raises(integrity.IntegrityError, match="Staged bytes/mode"):
        integrity.write_github_payload(repo, tmp_path / "filtered", source, source["head"])
    assert not (tmp_path / "filtered/candidate-metadata.json").exists()


def test_no_change_payload_has_empty_verified_patch(repo, tmp_path):
    source = integrity.snapshot(repo)
    result = integrity.write_github_payload(repo, tmp_path / "unchanged", source, source["head"])
    assert result["candidate_git_tree"] == git(repo, "rev-parse", "HEAD^{tree}").decode().strip()
    assert (tmp_path / "unchanged/candidate.patch").read_bytes() == b""
