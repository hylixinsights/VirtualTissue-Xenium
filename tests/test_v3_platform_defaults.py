"""Platform defaults must affect v3 while preserving explicit geometry."""
from virtualtissue_xenium.v2.context import validate_context
from virtualtissue_xenium.v2 import dynamics
from virtualtissue_xenium.v3 import ontology


def context(**overrides):
    value = dict(schema_version="3.0", tissue="skin", organism="human", disease="control",
                 assay="xenium", sampling=dict(design="subsampled", description="test", missing_populations=[]),
                 perturbations=[], recognition=False, time_step_hours=.5)
    value.update(overrides)
    return value


def changed_profile(monkeypatch):
    original = ontology.load_platform

    def load(name):
        profile = original(name)
        profile.update(neighbor_radius_um=73.0, mean_radius_um=9.0)
        return profile

    monkeypatch.setattr(ontology, "load_platform", load)


def test_v3_geometric_defaults_come_from_platform_with_context_precedence(monkeypatch):
    changed_profile(monkeypatch)
    declared = validate_context(context())
    assert declared["neighbor_radius_um"] == 73
    assert declared["dynamics"]["default_radius_um"] == 9
    explicit = validate_context(context(neighbor_radius_um=41, dynamics={"default_radius_um": 7}))
    assert explicit["neighbor_radius_um"] == 41
    assert explicit["dynamics"]["default_radius_um"] == 7


def test_lineage_radius_precedes_platform_fallback_and_is_not_bin_size(monkeypatch):
    changed_profile(monkeypatch)
    declared = validate_context(context())
    cell = {"lineage": "T_CD8"}
    assert dynamics.radius(cell, declared) == ontology.lineage_rule("T_CD8")["mean_radius_um"]
    monkeypatch.setattr(dynamics, "lineage_rule", lambda cell, context: {})
    assert dynamics.radius(cell, declared) == 9


def test_v2_geometric_default_does_not_load_v3_platform(monkeypatch):
    def unexpected(name):
        raise AssertionError("Schema 2.0 must not consult v3 platform defaults")

    monkeypatch.setattr(ontology, "load_platform", unexpected)
    declared = validate_context(context(schema_version="2.0"))
    assert declared["neighbor_radius_um"] == 30
    assert "dynamics" not in declared
