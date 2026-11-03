"""Native templates use measured geometry, never expression or inferred labels."""
import csv
import json

import h5py
import numpy as np
import pandas as pd
import pytest

from virtualtissue_xenium.v3.adapters import write_template
from virtualtissue_xenium.v3.adapters.common import annotation_rows, bounded_coordinate_records


def rows(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle))


def assert_blank(records, unit):
    assert records
    assert {r['measurement_unit'] for r in records} == {unit}
    assert all(r['cell_type'] == r['lineage'] == r['purity'] == '' for r in records)


def test_atera_template_needs_only_obs_geometry_and_is_not_an_annotation(tmp_path):
    import zarr
    from zarr.storage import ZipStore
    with ZipStore(str(tmp_path/'cell_feature_matrix.zarr.zip'), mode='w') as store:
        v3 = int(zarr.__version__.split('.')[0]) >= 3
        root = zarr.open_group(store=store, mode='w', **({'zarr_format': 2} if v3 else {}))
        for key, values in {'obs/barcode': ['b', 'a'], 'obs/centroid_column': [20., 10.],
                            'obs/centroid_row': [4., 2.]}.items():
            (root.create_array if v3 else root.create_dataset)(key, data=np.asarray(values))
        # X and var are deliberately absent: templates must never touch them.
    output = write_template('atera', tmp_path, tmp_path/'review.csv')
    annotations = rows(output)
    assert_blank(annotations, 'segmented_cell')
    assert [(r['cell_id'], float(r['x']), float(r['y'])) for r in annotations] == [('a', 10., 2.), ('b', 20., 4.)]
    with pytest.raises(ValueError, match='explicitly name'):
        annotation_rows(output, {'a', 'b'})
    saved = output.read_bytes()
    with pytest.raises(FileExistsError):
        write_template('atera', tmp_path, output)
    assert output.read_bytes() == saved


def test_visium_template_reads_barcodes_and_scales_without_matrix_values(tmp_path):
    spatial = tmp_path/'spatial'
    spatial.mkdir()
    (spatial/'scalefactors_json.json').write_text(json.dumps(dict(microns_per_pixel=.5, bin_size_um=8)))
    pd.DataFrame(dict(barcode=['a', 'b', 'c'], pxl_col_in_fullres=[0., 16., 32.],
                      pxl_row_in_fullres=[0., 0., 0.])).to_parquet(spatial/'tissue_positions.parquet')
    with h5py.File(tmp_path/'filtered_feature_bc_matrix.h5', 'w') as handle:
        handle.create_dataset('matrix/barcodes', data=np.array(['c', 'a', 'b'], dtype='S'))
        # No counts, indices, pointers or features are present.
    output = write_template('visium_hd', tmp_path, tmp_path/'review.csv',
                            roi=dict(center=[8., 0.], radius_um=2., max_cells=1))
    annotations = rows(output)
    assert_blank(annotations, 'binned_spot')
    assert [(r['cell_id'], float(r['x'])) for r in annotations] == [('b', 8.)]
    with pytest.raises(ValueError, match='explicitly name'):
        annotation_rows(output, {'a', 'b', 'c'})
    with pytest.raises(ValueError, match='microns_per_pixel'):
        write_template('visium_hd', tmp_path, tmp_path/'invalid.csv', measurement={'microns_per_pixel': float('nan')})
    assert not (tmp_path/'invalid.csv').exists()


def test_gem_bin_template_streams_unique_positions_and_keeps_nearest_ties(tmp_path):
    gem = tmp_path/'bins.gem'
    # Geometry is sufficient; no expression column is read or accumulated.
    gem.write_text('#BinType=Bin\nx\ty\n100\t0\n2\t0\n1\t0\n0\t0\n1\t0\n-1\t0\n')
    measurement = dict(unit='binned_spot', microns_per_coordinate=.5)
    output = write_template('stereo_seq', gem, tmp_path/'review.csv',
                            roi=dict(center=[0., 0.], max_cells=2), measurement=measurement)
    annotations = rows(output)
    assert_blank(annotations, 'binned_spot')
    assert [(r['cell_id'], float(r['x'])) for r in annotations] == [('-1:0', -.5), ('0:0', 0.)]
    with pytest.raises(ValueError, match='explicit center'):
        write_template('stereo_seq', gem, tmp_path/'unbounded.csv',
                       roi=dict(max_cells=2), measurement=measurement)
    assert not (tmp_path/'unbounded.csv').exists()


def test_streaming_nearest_selection_matches_complete_geometry_selection():
    # Repeats and late nearer positions must not consume the bounded candidate slots.
    from virtualtissue_xenium.v3.adapters.common import select_indices
    positions = [(f'c{i:04d}', float(i), float(i % 3)) for i in range(200)]
    streamed = bounded_coordinate_records(iter(positions[::-1] + positions),
                                          dict(center=[31., 0.], max_cells=7, radius_um=15.))
    selected = select_indices([r[0] for r in positions], [(r[1], r[2]) for r in positions],
                              center=[31., 0.], max_cells=7, radius_um=15.)
    assert [r['cell_id'] for r in streamed] == [positions[i][0] for i in selected]


def test_segmented_gem_requires_external_measured_centroids_and_exact_ids(tmp_path):
    gem = tmp_path/'cells.gem'
    gem.write_text('geneID\tx\ty\tMIDCount\tCellID\nMki67\t1000\t2000\t99\ta\nActb\t9000\t8000\t1\tb\n')
    measurement = dict(unit='segmented_cell', microns_per_coordinate=.5)
    with pytest.raises(ValueError, match='centroids_csv'):
        write_template('stereo_seq', gem, tmp_path/'invented.csv', measurement=measurement)
    centroids = tmp_path/'centroids.csv'
    centroids.write_text('cell_id,x,y\nb,30,40\na,10,20\n')
    measurement['centroids_csv'] = str(centroids)
    output = write_template('stereo_seq', gem, tmp_path/'review.csv', measurement=measurement)
    annotations = rows(output)
    assert_blank(annotations, 'segmented_cell')
    assert [(r['cell_id'], float(r['x']), float(r['y'])) for r in annotations] == [('a', 10., 20.), ('b', 30., 40.)]
    assert float(annotations[0]['x']) != 1000 * .5  # Not the transcript coordinate.
    centroids.write_text('cell_id,x,y\nmissing,10,20\n')
    with pytest.raises(ValueError, match='absent from GEM'):
        write_template('stereo_seq', gem, tmp_path/'missing.csv', measurement=measurement)
    assert not (tmp_path/'missing.csv').exists()


def test_bad_coordinate_stream_fails_before_writing_template(tmp_path):
    gem = tmp_path/'bad.gem'
    gem.write_text('x\ty\nnan\t0\n')
    with pytest.raises(ValueError, match='finite'):
        write_template('stereo_seq', gem, tmp_path/'bad.csv',
                       measurement=dict(unit='binned_spot', microns_per_coordinate=.5))
    assert not (tmp_path/'bad.csv').exists()
