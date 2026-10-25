"""Download the public example sample from 10x Genomics.

The example is *Xenium Prime Human Lung Cancer FFPE, experiment 2* (Xenium Prime
5K Human Pan Tissue and Pathways panel, Xenium Onboard Analysis 3.0.0), published
by 10x Genomics under CC BY 4.0:
https://www.10xgenomics.com/datasets/xenium-human-lung-cancer-post-xenium-technote

10x distributes it as one 22.6 GB ZIP ("Xenium Output Bundle"). Part 1 needs only
ten files from that ZIP (about 503 MiB). :func:`download_public_lung5k` reads the
remote ZIP's directory and transfers just those ten members with HTTP byte-range
requests, then verifies each file against the SHA-256 values recorded when the
example was first prepared (``datasets/xenium_prime_lung5k/DOWNLOAD_MANIFEST.json``).

If you prefer, download the full ZIP yourself, extract the ten files into
``<data_dir>/raw/`` and call :func:`register_manual_download`.

Only the Python standard library is used here.
"""
from __future__ import annotations

import hashlib
import io
import json
import shutil
import time
import urllib.error
import urllib.request
import zipfile
from pathlib import Path

from . import DATASETS_DIR

PUBLIC_SAMPLE = "xenium_prime_lung5k"
PUBLIC_SAMPLE_ID = "Xenium_Prime_Human_Lung_Cancer_FFPE"


def dataset_metadata(name: str = PUBLIC_SAMPLE) -> tuple[dict, dict]:
    """Return ``(SOURCE.json, DOWNLOAD_MANIFEST.json)`` of a packaged example dataset."""
    folder = DATASETS_DIR / name
    source = json.loads((folder / "SOURCE.json").read_text())
    manifest = json.loads((folder / "DOWNLOAD_MANIFEST.json").read_text())
    return source, manifest


def sha256_file(path, chunk: int = 8 * 1024 * 1024) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


class RemoteZipFile(io.RawIOBase):
    """Read-only, seekable view of a remote file through HTTP ``Range`` requests."""

    def __init__(self, url: str, *, timeout: float = 120, retries: int = 5, pause_s: float = 3.0):
        super().__init__()
        self.url, self.timeout, self.retries, self.pause_s = url, timeout, retries, pause_s
        self.position = 0
        self.requests = 0
        self.size = self._probe_size()

    def _open(self, start: int, end: int):
        request = urllib.request.Request(self.url, headers={"Range": f"bytes={start}-{end}",
                                                            "User-Agent": "virtualtissue-xenium"})
        return urllib.request.urlopen(request, timeout=self.timeout)

    def _probe_size(self) -> int:
        with self._open(0, 0) as response:
            content_range = response.headers.get("Content-Range", "")
            if response.status != 206 or "/" not in content_range:
                raise RuntimeError("The server does not support byte-range requests; "
                                   "download the full ZIP and use register_manual_download().")
            return int(content_range.rsplit("/", 1)[1])

    def _get(self, start: int, end: int) -> bytes:
        last_error = None
        for attempt in range(1, self.retries + 1):
            try:
                with self._open(start, end) as response:
                    if response.status != 206:
                        raise RuntimeError(f"Expected HTTP 206, received {response.status}.")
                    data = response.read()
                self.requests += 1
                if len(data) != end - start + 1:
                    raise RuntimeError("Incomplete byte range received.")
                return data
            except (urllib.error.URLError, TimeoutError, ConnectionError, RuntimeError, OSError) as error:
                last_error = error
                if attempt < self.retries:
                    time.sleep(self.pause_s * attempt)
        raise RuntimeError(f"Download failed after {self.retries} attempts: {last_error}")

    # io.RawIOBase interface -------------------------------------------------
    def readable(self) -> bool:
        return True

    def seekable(self) -> bool:
        return True

    def tell(self) -> int:
        return self.position

    def seek(self, offset: int, whence: int = io.SEEK_SET) -> int:
        if whence == io.SEEK_SET:
            position = offset
        elif whence == io.SEEK_CUR:
            position = self.position + offset
        elif whence == io.SEEK_END:
            position = self.size + offset
        else:
            raise ValueError("Invalid whence.")
        if position < 0:
            raise ValueError("Negative seek position.")
        self.position = position
        return position

    def readinto(self, buffer) -> int:
        if self.position >= self.size:
            return 0
        n = min(len(buffer), self.size - self.position)
        data = self._get(self.position, self.position + n - 1)
        buffer[:n] = data
        self.position += n
        return n


def _copy_verified(source, target: Path, expected_sha256: str, expected_bytes: int, progress=print) -> None:
    temporary = target.with_name(target.name + ".partial")
    digest, written, next_report = hashlib.sha256(), 0, 50 * 1024 * 1024
    with temporary.open("wb") as out:
        while True:
            block = source.read(4 * 1024 * 1024)
            if not block:
                break
            out.write(block)
            digest.update(block)
            written += len(block)
            if progress and written >= next_report:
                progress(f"    {written / 2**20:,.0f} / {expected_bytes / 2**20:,.0f} MiB")
                next_report += 50 * 1024 * 1024
    if written != expected_bytes or digest.hexdigest() != expected_sha256:
        temporary.unlink(missing_ok=True)
        raise RuntimeError(f"Checksum mismatch for {target.name}; the partial file was removed.")
    temporary.replace(target)


def download_from_archive(archive_url: str, manifest: dict, source: dict, data_dir, *,
                          block_size: int = 8 * 1024 * 1024, progress=print) -> Path:
    """Extract the members listed in ``manifest`` from a remote ZIP into ``data_dir``.

    Files that already exist with the expected SHA-256 are kept. ``SOURCE.json`` and
    ``DOWNLOAD_MANIFEST.json`` are written next to ``raw/`` at the end, so a partially
    downloaded folder is never mistaken for a complete one.
    """
    data_dir = Path(data_dir)
    (data_dir / "raw").mkdir(parents=True, exist_ok=True)
    pending = []
    for item in manifest["files"]:
        target = data_dir / item["path"]
        if target.is_file() and target.stat().st_size == item["bytes"] and sha256_file(target) == item["sha256"]:
            progress(f"  already present and verified: {item['path']}")
        else:
            pending.append(item)
    if pending:
        total = sum(i["bytes"] for i in pending)
        progress(f"Downloading {len(pending)} file(s), {total / 2**20:,.1f} MiB, from {archive_url}")
        remote = RemoteZipFile(archive_url)
        with io.BufferedReader(remote, buffer_size=block_size) as reader, zipfile.ZipFile(reader) as archive:
            for item in pending:
                progress(f"  {item['path']} ({item['bytes'] / 2**20:,.1f} MiB)")
                info = archive.getinfo(item["source_archive_member"])
                with archive.open(info) as member:  # zipfile also checks the member CRC32
                    _copy_verified(member, data_dir / item["path"], item["sha256"], item["bytes"], progress)
        progress(f"  HTTP range requests: {remote.requests}")
    (data_dir / "SOURCE.json").write_text(json.dumps(source, indent=2) + "\n")
    (data_dir / "DOWNLOAD_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")
    progress(f"Complete and verified: {data_dir}")
    return data_dir


def download_public_lung5k(data_dir, *, progress=print) -> Path:
    """Download the ten original files of the public example into ``data_dir``.

    Resulting layout (exactly what Part 1 expects as ``DATA_DIR``)::

        data_dir/
            SOURCE.json
            DOWNLOAD_MANIFEST.json
            raw/  (10 unmodified 10x files, ~503 MiB)
    """
    source, manifest = dataset_metadata(PUBLIC_SAMPLE)
    return download_from_archive(source["archive_url"], manifest, source, data_dir, progress=progress)


def verify_data_folder(data_dir, *, check_hashes: bool = True) -> dict:
    """Check that every file listed in ``DOWNLOAD_MANIFEST.json`` is present and unchanged."""
    data_dir = Path(data_dir)
    manifest = json.loads((data_dir / "DOWNLOAD_MANIFEST.json").read_text())
    problems = []
    for item in manifest["files"]:
        path = data_dir / item["path"]
        if not path.is_file():
            problems.append(f"missing: {item['path']}")
        elif path.stat().st_size != item["bytes"]:
            problems.append(f"size differs: {item['path']}")
        elif check_hashes and sha256_file(path) != item["sha256"]:
            problems.append(f"checksum differs: {item['path']}")
    return dict(data_dir=str(data_dir), files=len(manifest["files"]), complete=bool(manifest.get("complete")),
                hashes_checked=check_hashes, problems=problems, ok=not problems and bool(manifest.get("complete")))


def register_manual_download(data_dir, *, progress=print) -> dict:
    """Use files you extracted yourself from the full 10x ZIP.

    Put the ten files (see ``DOWNLOAD_MANIFEST.json``) into ``<data_dir>/raw/``; this
    writes the manifest and source description and verifies every checksum.
    """
    data_dir = Path(data_dir)
    source, manifest = dataset_metadata(PUBLIC_SAMPLE)
    for name, obj in (("SOURCE.json", source), ("DOWNLOAD_MANIFEST.json", manifest)):
        (data_dir / name).write_text(json.dumps(obj, indent=2) + "\n")
    report = verify_data_folder(data_dir)
    if report["problems"]:
        progress("Problems found:\n  " + "\n  ".join(report["problems"]))
    else:
        progress(f"All {report['files']} files verified in {data_dir}")
    return report


def copy_packaged_checksums(data_dir) -> Path:
    """Write ``CHECKSUMS.sha256`` (for ``sha256sum -c``) into ``data_dir``."""
    target = Path(data_dir) / "CHECKSUMS.sha256"
    shutil.copyfile(DATASETS_DIR / PUBLIC_SAMPLE / "CHECKSUMS.sha256", target)
    return target
