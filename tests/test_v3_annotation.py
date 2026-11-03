"""ROI rules and annotation agreement must not manufacture reviewed labels."""
import csv

import pytest

from virtualtissue_xenium.v3.annotation import compare_reviewers, select_rois, write_annotation_comparison


def tissue(sample="s"):
    labels = ["T_CD8", "malignant", "malignant", "malignant"] + ["malignant"] * 4
    return [dict(cell_id=f"{sample}-{i}", sample_id=sample, lineage=label, cell_type=label,
                 x=float((0 if i < 4 else 1000) + i % 4), y=0.) for i, label in enumerate(labels)]


def write_csv(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["cell_id", "lineage"])
        writer.writeheader()
        writer.writerows(dict(cell_id=cid, lineage=lineage) for cid, lineage in rows)
    return path


def test_roi_rule_finds_front_core_and_is_order_invariant():
    cells = tissue()
    first = select_rois(cells, cells_per_roi=4)
    assert first == select_rois(list(reversed(cells)), cells_per_roi=4)
    assert first["complete"]
    result = first["samples"]["s"]
    assert result["invasive_front"]["center_um"] == [0., 0.]
    assert result["invasive_front"]["cd8_malignant_contacts"] == 3
    assert result["tumour_core"]["center_um"] == [1000., 0.]
    assert result["tumour_core"]["immune_fraction"] == 0
    assert not set(result["invasive_front"]["cell_ids"]) & set(result["tumour_core"]["cell_ids"])


def test_samples_never_share_rois_and_candidate_limit_is_deterministic():
    result = select_rois(tissue("a") + tissue("b"), cells_per_roi=4, max_candidates=1)
    assert set(result["samples"]) == {"a", "b"}
    for sample, selections in result["samples"].items():
        assert len(selections["candidates"]) == 1
        assert all(cid.startswith(sample) for cid in selections["invasive_front"]["cell_ids"])
        assert selections["tumour_core"] is None
        assert not selections["complete"]


def test_unknown_or_missing_geometry_cannot_be_an_annotated_roi():
    cells = tissue()
    for c in cells:
        c["lineage"] = "unknown"
    with pytest.raises(ValueError, match="reviewed malignant and T_CD8"):
        select_rois(cells, cells_per_roi=4)
    cells = tissue()
    del cells[0]["x"]
    with pytest.raises(ValueError, match="measured x/y"):
        select_rois(cells, cells_per_roi=4)
    with pytest.raises(ValueError, match="smaller ROI"):
        select_rois(tissue(), cells_per_roi=1000)


def test_unavailable_roi_is_null_and_threshold_is_not_relaxed():
    cells = tissue()[:4]
    result = select_rois(cells, cells_per_roi=4, require_disjoint=False, core_max_immune_fraction=.25)
    assert result["samples"]["s"]["tumour_core"] is None  # exactly 25%, but rule is strictly below
    assert not result["complete"]
    result = select_rois(tissue(), cells_per_roi=4, contact_radius_um=.01)
    assert result["samples"]["s"]["invasive_front"] is None
    assert "no_CD8_malignant_contacts_at_declared_radius" in result["samples"]["s"]["missing_reasons"]


def test_kappa_known_table_and_export(tmp_path):
    a = write_csv(tmp_path / "a.csv", [("1", "A"), ("2", "A"), ("3", "B"), ("4", "B")])
    b = write_csv(tmp_path / "b.csv", [("4", "B"), ("2", "B"), ("1", "A"), ("3", "B")])
    report = write_annotation_comparison(a, b, tmp_path / "agreement")
    assert report["overall"]["observed_agreement"] == .75
    assert report["overall"]["kappa"] == .5
    assert report["by_lineage"]["A"]["kappa"] == .5
    assert report["disagreements"] == [dict(cell_id="2", reviewer_1="A", reviewer_2="B")]
    assert (tmp_path / "agreement/disagreements.csv").read_text().splitlines()[1] == "2,A,B"


def test_mismatched_ids_require_explicit_intersection(tmp_path):
    a = write_csv(tmp_path / "a.csv", [("1", "A"), ("2", "B")])
    b = write_csv(tmp_path / "b.csv", [("1", "A"), ("3", "B")])
    with pytest.raises(ValueError, match="IDs differ"):
        compare_reviewers(a, b)
    report = compare_reviewers(a, b, id_scope="intersection")
    assert report["common_cell_ids"] == ["1"]
    assert report["excluded_only_first"] == ["2"] and report["excluded_only_second"] == ["3"]
    assert report["overall"]["kappa"] is None
    assert report["overall"]["observed_agreement"] == 1


def test_no_common_ids_and_constant_agreement_have_null_kappa(tmp_path):
    a = write_csv(tmp_path / "a.csv", [("1", "unknown")])
    b = write_csv(tmp_path / "b.csv", [("2", "unknown")])
    report = compare_reviewers(a, b, id_scope="intersection")
    assert report["overall"]["kappa"] is None
    assert report["overall"]["kappa_status"] == "undefined_no_common_cells"
    report = compare_reviewers(a, a)
    assert report["overall"]["kappa"] is None
    assert report["overall"]["kappa_status"] == "undefined_zero_expected_disagreement"


def test_duplicate_or_blank_reviewer_ids_are_rejected(tmp_path):
    a = write_csv(tmp_path / "a.csv", [("1", "A"), ("1", "B")])
    with pytest.raises(ValueError, match="Duplicate"):
        compare_reviewers(a, a)
    a = write_csv(tmp_path / "a.csv", [("1", "")])
    with pytest.raises(ValueError, match="must be filled"):
        compare_reviewers(a, a)
