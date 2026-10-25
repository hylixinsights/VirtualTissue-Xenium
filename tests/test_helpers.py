"""Tests of the new helper modules (download, own data, Drive backup, results, CLI)."""
import hashlib
import http.server
import io
import json
import os
import shutil
import tarfile
import tempfile
import threading
import unittest
import zipfile
from pathlib import Path

import virtualtissue_xenium as vtx
from virtualtissue_xenium import colab, download, own_data

BUNDLE = os.environ.get("VTX_TEST_BUNDLE")  # optional: path to a simulation_input_*.zip


class RangeHandler(http.server.BaseHTTPRequestHandler):
    payload = b""

    def log_message(self, *args):
        pass

    def do_GET(self):
        data = self.payload
        header = self.headers.get("Range")
        if not header:
            self.send_response(200)
            self.send_header("Content-Length", str(len(data)))
            self.end_headers()
            self.wfile.write(data)
            return
        start, end = header.split("=", 1)[1].split("-")
        start, end = int(start), min(int(end), len(data) - 1)
        chunk = data[start:end + 1]
        self.send_response(206)
        self.send_header("Content-Range", f"bytes {start}-{end}/{len(data)}")
        self.send_header("Content-Length", str(len(chunk)))
        self.end_headers()
        self.wfile.write(chunk)


class DownloadTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        members = {"a.bin": os.urandom(3_000_000), "b.txt": b"hello xenium\n" * 1000, "unused.bin": os.urandom(10_000)}
        buffer = io.BytesIO()
        with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as z:
            for name, content in members.items():
                z.writestr(name, content)
        RangeHandler.payload = buffer.getvalue()
        self.server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), RangeHandler)
        threading.Thread(target=self.server.serve_forever, daemon=True).start()
        self.url = f"http://127.0.0.1:{self.server.server_address[1]}/archive.zip"
        self.manifest = dict(sample_id="test", complete=True, files=[
            dict(path=f"raw/{n}", bytes=len(members[n]), sha256=hashlib.sha256(members[n]).hexdigest(),
                 source_archive_member=n) for n in ("a.bin", "b.txt")])
        self.source = dict(sample_id="test", archive_url=self.url)

    def tearDown(self):
        self.server.shutdown()
        self.temp.cleanup()

    def test_range_download_verifies_and_resumes(self):
        dest = Path(self.temp.name) / "data"
        download.download_from_archive(self.url, self.manifest, self.source, dest, block_size=256 * 1024, progress=lambda m: None)
        self.assertTrue((dest / "DOWNLOAD_MANIFEST.json").is_file())
        self.assertFalse((dest / "raw" / "unused.bin").exists())
        self.assertTrue(download.verify_data_folder(dest)["ok"])
        messages = []
        download.download_from_archive(self.url, self.manifest, self.source, dest, progress=messages.append)
        self.assertTrue(all("already present" in m for m in messages[:2]))

    def test_checksum_mismatch_is_rejected(self):
        bad = json.loads(json.dumps(self.manifest))
        bad["files"][0]["sha256"] = "0" * 64
        with self.assertRaisesRegex(RuntimeError, "Checksum mismatch"):
            download.download_from_archive(self.url, bad, self.source, Path(self.temp.name) / "bad", progress=lambda m: None)
        self.assertFalse((Path(self.temp.name) / "bad" / "DOWNLOAD_MANIFEST.json").exists())

    def test_packaged_public_metadata(self):
        source, manifest = download.dataset_metadata()
        self.assertEqual(source["sample_id"], download.PUBLIC_SAMPLE_ID)
        self.assertEqual(len(manifest["files"]), 10)
        self.assertTrue(manifest["complete"])


def write_fake_xenium(raw: Path, n_cells=60, clusters=3):
    """Tiny but structurally valid Xenium output (10x HDF5, parquet tables, analysis.tar.gz)."""
    import h5py
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    raw.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(0)
    genes = [f"G{i}" for i in range(8)]
    ids = [f"cell{i:04d}-1" for i in range(n_cells)]
    X = sp.random(len(genes) + 2, n_cells, density=.5, format="csc", random_state=0, data_rvs=lambda n: rng.integers(1, 5, n)).astype(np.int32)
    with h5py.File(raw / "cell_feature_matrix.h5", "w") as f:
        m = f.create_group("matrix")
        m["barcodes"] = np.array(ids, dtype="S")
        m["data"], m["indices"], m["indptr"], m["shape"] = X.data, X.indices, X.indptr, np.array(X.shape)
        feats = m.create_group("features")
        feats["name"] = np.array(genes + ["NegCtrl1", "NegCtrl2"], dtype="S")
        feats["id"] = np.array(genes + ["N1", "N2"], dtype="S")
        feats["feature_type"] = np.array(["Gene Expression"] * len(genes) + ["Negative Control Probe"] * 2, dtype="S")
    pd.DataFrame(dict(cell_id=ids, x_centroid=rng.uniform(0, 100, n_cells), y_centroid=rng.uniform(0, 100, n_cells),
                      transcript_counts=np.asarray(X[:len(genes)].sum(axis=0)).ravel(), cell_area=50., nucleus_area=20.)).to_parquet(raw / "cells.parquet")
    square = [(0, 0), (1, 0), (1, 1), (0, 1), (0, 0)]
    poly = pd.DataFrame([dict(cell_id=c, vertex_x=x, vertex_y=y, label_id=i) for i, c in enumerate(ids) for x, y in square])
    poly.to_parquet(raw / "cell_boundaries.parquet")
    poly.to_parquet(raw / "nucleus_boundaries.parquet")
    cl = pd.DataFrame(dict(Barcode=ids, Cluster=[1 + i % clusters for i in range(n_cells)]))
    de = pd.DataFrame({"Feature ID": genes, "Feature Name": genes})
    for c in range(1, clusters + 1):
        de[f"Cluster {c} Mean Counts"] = rng.uniform(0, 1, len(genes))
        de[f"Cluster {c} Log2 fold change"] = rng.normal(size=len(genes))
        de[f"Cluster {c} Adjusted p value"] = .001
    with tarfile.open(raw / "analysis.tar.gz", "w:gz") as t:
        for name, frame in ((own_data.CLUSTERS_MEMBER, cl), (own_data.DIFFEXP_MEMBER, de)):
            data = frame.to_csv(index=False).encode()
            info = tarfile.TarInfo(name)
            info.size = len(data)
            t.addfile(info, io.BytesIO(data))
    return ids


class OwnDataTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.outs = Path(self.temp.name) / "outs"
        write_fake_xenium(self.outs)

    def tearDown(self):
        self.temp.cleanup()

    def test_create_folder_template_validate_and_activate(self):
        data = Path(self.temp.name) / "data"
        own_data.create_data_folder(self.outs, data, "MY_SAMPLE", progress=lambda m: None)
        manifest = json.loads((data / "DOWNLOAD_MANIFEST.json").read_text())
        self.assertTrue(manifest["complete"])
        self.assertEqual({f["path"] for f in manifest["files"]}, {f"raw/{n}" for n in own_data.REQUIRED_FILES})
        self.assertTrue(download.verify_data_folder(data)["ok"])
        with self.assertRaisesRegex(ValueError, "not the public example"):
            own_data.require_annotation_map(data, "")
        template = own_data.write_annotation_template(data, data / "clusters.json")
        records = json.loads(template.read_text())
        self.assertEqual([r["cluster"] for r in records], [1, 2, 3])
        self.assertTrue(any("placeholder" in p for p in own_data.validate_annotations(template, data)))
        for r, (label, comp, runtime) in zip(records, [("Tumour", "tumor_candidate", "Tumor"), ("CAF", "fibroblast", "CAF"),
                                                       ("CD8 T", "immune", "T_CD8_cytotoxic")]):
            r.update(label=label, compartment=comp, runtime_label=runtime, markers=["G1"])
        template.write_text(json.dumps(records))
        self.assertEqual(own_data.validate_annotations(template, data), [])
        prep = vtx.prepare
        try:
            target = own_data.use_cluster_annotations(prep, template, Path(self.temp.name) / "work")
            self.assertEqual(prep.RESOURCES, target)
            self.assertEqual(json.loads((target / "cluster_annotations.json").read_text()), records)
            self.assertTrue((target / "go20_msigdb_2026.1.Hs.json").is_file())
            # The engine itself reads the custom map: load + annotate the fake sample.
            sample = prep.load_sample(data, Path(self.temp.name) / "out", verify_hashes=True, min_counts=1, min_genes=1)
            table = prep.annotate(sample, {})
            self.assertEqual(sorted(table.label), sorted(["Tumour", "CAF", "CD8 T"]))
        finally:
            own_data.use_packaged_annotations(prep)
        self.assertEqual(prep.RESOURCES, vtx.ENGINE_DIR / "resources")

    def test_missing_file_is_explained(self):
        (self.outs / "analysis.tar.gz").unlink()
        with self.assertRaisesRegex(ValueError, "tar -czf"):
            own_data.check_xenium_files(self.outs)


class BackupTests(unittest.TestCase):
    def test_backup_restore_skips_lock_and_temporary_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            run = Path(tmp) / "session" / "R03_run"
            run.mkdir(parents=True)
            for name in ("run_config.json", "jev_journal.jsonl", "episode.json.gz", ".running.lock", "x.tmp"):
                (run / name).write_text(name)
            (run.parent / "R03_run_results.zip").write_text("zip")
            drive = Path(tmp) / "drive"
            messages = []
            progress = colab.progress_with_backup(run, drive, printer=messages.append)
            progress("Jev requests recorded: 10/600")
            self.assertFalse(drive.exists())
            progress("Saved 1.00 simulated hours; 50 cell decisions.")
            self.assertEqual(sorted(p.name for p in (drive / "R03_run").iterdir()),
                             ["episode.json.gz", "jev_journal.jsonl", "run_config.json"])
            self.assertTrue((drive / "R03_run_results.zip").is_file())
            fresh = Path(tmp) / "new_session" / "R03_run"
            self.assertTrue(colab.restore_run(drive, fresh))
            self.assertEqual((fresh / "jev_journal.jsonl").read_text(), "jev_journal.jsonl")
            self.assertFalse(colab.restore_run(drive, fresh), "Never overwrite an existing local run.")

    def test_key_lookup_order_and_no_prompt(self):
        old = {k: os.environ.pop(k, None) for k in colab.KEY_NAMES}
        try:
            self.assertEqual(colab.get_jev_key(prompt=False), "")
            os.environ["JEV_API_KEY"] = "  second  "
            self.assertEqual(colab.get_jev_key(prompt=False), "second")
            os.environ["TYPESAFE_API_KEY"] = "first"
            self.assertEqual(colab.get_jev_key(prompt=False), "first")
        finally:
            for k, v in old.items():
                os.environ.pop(k, None)
                if v is not None:
                    os.environ[k] = v


@unittest.skipUnless(BUNDLE, "set VTX_TEST_BUNDLE=/path/to/simulation_input_*.zip")
class BundleRunTests(unittest.TestCase):
    def test_cli_check_run_local_and_explore(self):
        from virtualtissue_xenium import cli, results
        with tempfile.TemporaryDirectory() as tmp:
            cli.main(["check", BUNDLE, "--workdir", f"{tmp}/check", "--hours", "1"])
            self.assertTrue(Path(f"{tmp}/check/preflight/first_request_preview.json").is_file())
            cli.main(["run", BUNDLE, f"{tmp}/results", "--run-name", "cli_test", "--mode", "local", "--hours", "1"])
            zip_path = Path(f"{tmp}/results/cli_test_local_preview_results.zip")
            self.assertTrue(zip_path.is_file())
            res = results.load_results(zip_path)
            s = results.summary(res)
            self.assertTrue(s["completed"])
            self.assertEqual(s["simulated_hours"], 1.0)
            self.assertGreater(s["cell_decisions"], 0)
            self.assertEqual(int(results.population_by_type(res).iloc[0].sum()), res["episode"]["summary"]["final"]["acting"]["alive_in_roi"])
            self.assertFalse(results.choice_counts(res).empty)
            cli.main(["explore", str(zip_path), "--plot", f"{tmp}/overview.png"])
            self.assertTrue(Path(f"{tmp}/overview.png").is_file())


if __name__ == "__main__":
    unittest.main()
