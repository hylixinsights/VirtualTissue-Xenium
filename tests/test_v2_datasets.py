"""Input failures and provenance checks for the portable v2 dataset path."""
import hashlib
import io
import json
import tarfile

import h5py
import numpy as np
import pandas as pd
import pytest

from virtualtissue_xenium.v2 import datasets as D


@pytest.fixture
def xenium(tmp_path):
    raw = tmp_path / 'raw'; raw.mkdir()
    pd.DataFrame({'cell_id': ['b', 'a', 'c'], 'x_centroid': [10., 0., 20.],
                  'y_centroid': [0., 0., 0.]}).to_csv(raw / 'cells.csv', index=False)
    with h5py.File(raw / 'cell_feature_matrix.h5', 'w') as f:
        m = f.create_group('matrix'); ft = m.create_group('features')
        m['barcodes'] = np.array(['a', 'b', 'c'], dtype='S')
        ft['name'] = np.array(['MKI67', 'Control', 'CD3D'], dtype='S')
        ft['feature_type'] = np.array(['Gene Expression', 'Negative Control Probe', 'Gene Expression'], dtype='S')
        # a: MKI67=2, Control=9; b: CD3D=5; c: MKI67=3, CD3D=4
        m['shape'] = [3, 3]; m['indptr'] = [0, 2, 3, 5]
        m['indices'] = [0, 1, 2, 0, 2]; m['data'] = [2., 9., 5., 3., 4.]
    annotations = tmp_path / 'annotations.csv'
    pd.DataFrame({'cell_id': ['c', 'a', 'b'], 'cell_type': ['KC', 'KC', 'T cell'],
                  'lineage': ['epithelial', 'epithelial', 'T_CD8']}).to_csv(annotations, index=False)
    return tmp_path, annotations


def test_exact_id_alignment_controls_and_counts(xenium):
    path, annotations = xenium
    cells, genes = D.load_xenium_cells(path, annotations)
    assert genes == ['MKI67', 'CD3D']
    assert [x['cell_id'] for x in cells] == ['a', 'b', 'c']
    assert cells[0]['x'] == 0. and cells[1]['x'] == 10.
    assert cells[0]['expression'] == {'MKI67': 2.}
    assert cells[1]['expression'] == {'CD3D': 5.}
    assert cells[1]['lineage'] == 'T_CD8'


def test_large_input_requires_explicit_roi(xenium):
    path, annotations = xenium
    with pytest.raises(ValueError, match='explicit center'):
        D.load_xenium_cells(path, annotations, max_cells=2)
    cells, _ = D.load_xenium_cells(path, annotations, max_cells=2, center=(0., 0.))
    assert [x['cell_id'] for x in cells] == ['a', 'b']
    cells, _ = D.load_xenium_cells(path, annotations, radius_um=5., center=(0., 0.))
    assert [x['cell_id'] for x in cells] == ['a']


def test_blank_template_does_not_create_annotations(xenium, tmp_path):
    path, _ = xenium
    template = D.write_annotation_template(path, tmp_path / 'review.csv')
    table = pd.read_csv(template, keep_default_na=False)
    assert set(table.cell_id) == {'a', 'b', 'c'}
    assert (table.cell_type == '').all() and (table.lineage == '').all()
    with pytest.raises(ValueError, match='Fill every'):
        D.load_xenium_cells(path, template)


@pytest.mark.parametrize('corruption', ['missing', 'extra', 'duplicate'])
def test_annotation_join_rejects_bad_ids(xenium, corruption):
    path, annotations = xenium
    table = pd.read_csv(annotations)
    if corruption == 'missing': table = table.iloc[:2]
    elif corruption == 'extra': table.loc[0, 'cell_id'] = 'unknown_id'
    else: table.loc[0, 'cell_id'] = 'a'
    table.to_csv(annotations, index=False)
    with pytest.raises(ValueError): D.load_xenium_cells(path, annotations)


def test_missing_measured_coordinate_id_is_an_error(xenium):
    path, annotations = xenium
    table = pd.read_csv(path / 'raw/cells.csv'); table.loc[0, 'cell_id'] = 'wrong'
    table.to_csv(path / 'raw/cells.csv', index=False)
    with pytest.raises(ValueError, match='exactly the same'):
        D.load_xenium_cells(path, annotations)


def test_negative_count_is_an_error(xenium):
    path, annotations = xenium
    with h5py.File(path / 'raw/cell_feature_matrix.h5', 'r+') as f: f['matrix/data'][0] = -2.
    with pytest.raises(ValueError, match='nonnegative'): D.load_xenium_cells(path, annotations)


def test_download_checks_size_hash_and_cache(tmp_path, monkeypatch):
    content = b'real count data'
    manifest = {'id': 'test', 'files': [{'path': 'matrix.txt', 'url': 'https://example.test/data',
                  'bytes': len(content), 'sha256': hashlib.sha256(content).hexdigest()}]}
    monkeypatch.setattr(D, 'get_dataset', lambda name: manifest)
    calls = []
    def opener(request, **kwargs): calls.append(request); return io.BytesIO(content)
    monkeypatch.setattr(D.urllib.request, 'urlopen', opener)
    D.download_dataset('test', tmp_path)
    assert D.verify_dataset('test', tmp_path)['ok']
    assert json.loads((tmp_path / 'DOWNLOAD_RECEIPT.json').read_text())['ok']
    D.download_dataset('test', tmp_path)
    assert len(calls) == 1
    with pytest.raises(ValueError, match='exceed max_bytes'): D.download_dataset('test', tmp_path, max_bytes=1)
    (tmp_path / 'raw/matrix.txt').write_bytes(b'bad')
    assert not D.verify_dataset('test', tmp_path)['ok']


def test_oversized_response_never_becomes_a_download(tmp_path, monkeypatch):
    manifest = {'id': 'test', 'files': [{'path': 'matrix.txt', 'url': 'https://example.test/data',
                 'bytes': 1, 'sha256': hashlib.sha256(b'x').hexdigest()}]}
    monkeypatch.setattr(D, 'get_dataset', lambda name: manifest)
    monkeypatch.setattr(D.urllib.request, 'urlopen', lambda *a, **kw: io.BytesIO(b'xxxxxx'))
    with pytest.raises(ValueError, match='exceeds declared size'): D.download_dataset('test', tmp_path)
    assert not (tmp_path / 'raw/matrix.txt').exists()
    assert not (tmp_path / 'raw/matrix.txt.partial').exists()
    assert not (tmp_path / 'DOWNLOAD_RECEIPT.json').exists()


def test_vendor_qc_omissions_stay_visible_in_template(xenium, tmp_path):
    path, _ = xenium
    content = b'Barcode,Cluster\na,1\nb,2\n'
    with tarfile.open(path / 'raw/analysis.tar.gz', 'w:gz') as tf:
        info = tarfile.TarInfo('analysis/clustering/gene_expression_graphclust/clusters.csv')
        info.size = len(content); tf.addfile(info, io.BytesIO(content))
    out = D.write_annotation_template(path, tmp_path / 'review.csv')
    table = pd.read_csv(out, keep_default_na=False).set_index('cell_id')
    assert set(table.index) == {'a', 'b', 'c'}
    assert table.loc['c', 'cluster_id'] == ''


def test_failed_redownload_invalidates_stale_success_receipt(tmp_path, monkeypatch):
    manifest = {'id': 'test', 'files': [{'path': 'matrix.txt', 'url': 'https://example.test/data',
                 'bytes': 1, 'sha256': hashlib.sha256(b'x').hexdigest()}]}
    monkeypatch.setattr(D, 'get_dataset', lambda name: manifest)
    monkeypatch.setattr(D.urllib.request, 'urlopen', lambda *a, **kw: io.BytesIO(b'corruption'))
    (tmp_path / 'DOWNLOAD_RECEIPT.json').write_text('{"ok":true}')
    with pytest.raises(ValueError): D.download_dataset('test', tmp_path)
    assert not (tmp_path / 'DOWNLOAD_RECEIPT.json').exists()
