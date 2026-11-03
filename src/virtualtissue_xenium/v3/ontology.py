"""Validated, versioned biological declarations; never inferred from disease names.

Data are packaged with the code so installed wheels and source checkouts use the
same ontology. Public helpers return independent copies. A cited source does not
make an assumed numerical rate a measurement of the user's tissue.
"""
from __future__ import annotations

from copy import deepcopy
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path

DATA_ROOT = Path(__file__).with_name("data")
SUPPORTED_VERSION = "3.0"
_ACTIONS = {"REST", "PROLIFERATE", "APOPTOSIS", "ACTIVATE", "EXHAUST", "RECOVER",
            "SECRETE", "CLEAR_CORPSE", "KILL", "MOVE"}


def _read(relative):
    path = DATA_ROOT / relative
    return json.loads(path.read_text(encoding="utf-8"))


def file_hash(relative):
    path = (DATA_ROOT / relative).resolve()
    if not path.is_relative_to(DATA_ROOT.resolve()) or not path.is_file():
        raise ValueError("Ontology data path must identify a packaged file.")
    return hashlib.sha256(path.read_bytes()).hexdigest()


def validate_ontology(document, kinetics):
    """Reject incomplete permissions and rates before they can reach a guidebook."""
    if document.get("ontology_version") != SUPPORTED_VERSION:
        raise ValueError("Unsupported ontology_version; expected '3.0'.")
    rows = document.get("lineages")
    if not isinstance(rows, dict) or not rows or "unknown" not in rows:
        raise ValueError("Ontology requires lineages including unknown.")
    hazards = kinetics.get("hazards", {})
    if set(hazards) != set(rows):
        raise ValueError("Ontology lineage and kinetics rows must match exactly.")
    for name, row in rows.items():
        actions = row.get("allowed_actions")
        if (row.get("name") != name or not isinstance(actions, list) or not actions
                or len(actions) != len(set(actions)) or "REST" not in actions
                or set(actions) - _ACTIONS):
            raise ValueError(f"Invalid allowed_actions for lineage {name}.")
        parent = row.get("parent")
        if parent is not None and (parent not in rows or parent == name):
            raise ValueError(f"Invalid parent for lineage {name}.")
        visited = {name}
        while parent is not None:
            if parent in visited:
                raise ValueError("Ontology parents contain a cycle.")
            visited.add(parent)
            parent = rows[parent].get("parent")
        radius = row.get("mean_radius_um")
        if type(radius) not in (int, float) or not math.isfinite(radius) or radius <= 0:
            raise ValueError(f"Invalid radius for lineage {name}.")
        for flag in ("motile", "can_kill", "can_be_killed_target", "can_clear", "proliferation_requires_activation"):
            if type(row.get(flag)) is not bool:
                raise ValueError(f"Ontology {name}.{flag} must be a boolean.")
        for flag, action in (("motile", "MOVE"), ("can_kill", "KILL"), ("can_clear", "CLEAR_CORPSE")):
            if row[flag] != (action in actions):
                raise ValueError(f"Ontology permission {name}.{flag} must match {action}.")
        if set(hazards[name]) != set(actions) - {"REST"}:
            raise ValueError(f"Kinetics actions do not match allowed_actions for {name}.")
        for action, rate in hazards[name].items():
            value = rate.get("hazard_per_hour")
            if type(value) not in (int, float) or not math.isfinite(value) or value < 0:
                raise ValueError(f"Invalid hazard for {name}.{action}.")
            if rate.get("confidence") not in kinetics.get("confidence_levels", {}):
                raise ValueError(f"Missing confidence for {name}.{action}.")
            if not rate.get("source"):
                raise ValueError(f"Missing source/assumption for {name}.{action}.")
    if rows["unknown"]["allowed_actions"] != ["REST"]:
        raise ValueError("Unknown lineage must remain REST-only.")


@lru_cache(maxsize=1)
def _ontology(version="3.0"):
    if version != SUPPORTED_VERSION:
        raise ValueError("Unsupported ontology_version; expected '3.0'.")
    document = _read("ontology/lineages_v3.json")
    kinetics = _read("ontology/kinetics_v3.json")
    validate_ontology(document, kinetics)
    return dict(document, kinetics=kinetics, hashes={
        "lineages_sha256": file_hash("ontology/lineages_v3.json"),
        "kinetics_sha256": file_hash("ontology/kinetics_v3.json"),
        "data_version_sha256": file_hash("DATA_VERSION"),
        "dataset_registry_sha256": file_hash("datasets/manifests.json"),
    })


def load_ontology(version="3.0"):
    return deepcopy(_ontology(version))


def lineage_rules(version="3.0"):
    return deepcopy(_ontology(version)["lineages"])


def lineage_rule(name, version="3.0"):
    """Read one lineage without repeatedly copying the entire rate table."""
    try:
        return deepcopy(_ontology(version)["lineages"][name])
    except KeyError as exc:
        raise ValueError(f"Unsupported lineage: {name!r}.") from exc


def load_signatures(organism):
    aliases = {"human": "human", "homo sapiens": "human", "mouse": "mouse", "mus musculus": "mouse"}
    species = aliases.get(str(organism).strip().casefold())
    if species is None:
        raise ValueError("No packaged signatures for this organism; provide explicit reviewed gene_sets.")
    filename = f"signatures/{species}_v3.json"
    source = _read(filename)
    gene_sets = {}
    for name, record in source["sets"].items():
        genes = record.get("genes")
        if (not isinstance(genes, list) or not genes or len(genes) != len(set(genes))
                or any(not isinstance(g, str) or not g.strip() or g != g.strip() for g in genes)):
            raise ValueError(f"Invalid signature membership: {name}.")
        gene_sets[name] = list(genes)
    provenance = deepcopy(source)
    provenance["file_sha256"] = file_hash(filename)
    if species == "mouse":
        provenance["ortholog_table_sha256"] = file_hash("signatures/human_mouse_orthologs_v3.json")
    return gene_sets, provenance


def load_platform(platform):
    if platform not in {"xenium", "scrna", "snrna", "atera", "visium_hd", "stereo_seq"}:
        raise ValueError(f"Unsupported platform: {platform!r}.")
    result = _read(f"platforms/{platform}.json")
    result["sha256"] = file_hash(f"platforms/{platform}.json")
    return result


def load_profile(name):
    if name not in {"skin_melanoma", "brain_cortex"}:
        raise ValueError(f"Unsupported tissue profile: {name!r}.")
    result = _read(f"profiles/{name}.json")
    result["sha256"] = file_hash(f"profiles/{name}.json")
    return result


def kinetics_for(context):
    ontology = load_ontology(context.get("ontology_version", "3.0"))
    table = ontology["kinetics"]
    profile = load_profile(context["tissue_profile"]) if context.get("tissue_profile") else None
    if profile:
        aliases = {"human": "human", "homo sapiens": "human", "mouse": "mouse", "mus musculus": "mouse"}
        if aliases.get(str(context.get("organism", "")).casefold()) != profile["organism"]:
            raise ValueError("Tissue profile organism does not match the declared organism.")
        for lineage, actions in profile["hazard_overrides"].items():
            if lineage not in table["hazards"]:
                raise ValueError("Tissue profile references an unsupported lineage.")
            for action, row in actions.items():
                if action not in table["hazards"][lineage]:
                    raise ValueError("Tissue profile cannot license a new action.")
                table["hazards"][lineage][action] = deepcopy(row)
        validate_ontology(ontology, table)
    hazards = {name: {action: float(row["hazard_per_hour"]) for action, row in rows.items()}
               for name, rows in table["hazards"].items()}
    return table, hazards, ontology, profile
