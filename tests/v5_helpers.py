"""Small synthetic inputs for the v5 tests; engineering fixtures, not biology."""
from copy import deepcopy

from virtualtissue_xenium.v2.inputs import prepare
from virtualtissue_xenium.v5.evidence import POLICY_VERSION

PANEL = ["MKI67", "TOP2A", "PCNA", "MCM2", "MCM5", "CDK1", "ACTB", "GAPDH", "CD69", "IL2RA", "IFNG", "TNFRSF9", "FOS", "JUN"]


def context(**overrides):
    result = dict(schema_version="3.0", tissue="lung", organism="human", disease="test", assay="xenium",
                  sampling=dict(design="subsampled", description="test", missing_populations=[]),
                  perturbations=[], recognition=False, time_step_hours=.5,
                  scoring=dict(method="depth_conditioned_auc", depth_bins=1, min_reference_cells=5, min_transcripts=0,
                               auc_top_fraction=1.0),
                  dynamics=dict(crowding_enabled=False))
    result.update(overrides)
    return result


def cells(n=12, *, scores=True):
    """n tumour cells and n T cells on a line; measured scores increase with the index."""
    result = []
    for i in range(n):
        for j, (label, lineage) in enumerate((("tumour", "malignant"), ("tcell", "T_CD8"))):
            cell = dict(cell_id=f"{label}{i:02}", cell_type=label, lineage=lineage, sample_id="s",
                        x=float(i) * 14, y=float(j) * 10,
                        expression={"MKI67": i + 1., "TOP2A": 2., "PCNA": 3., "MCM2": 4., "MCM5": 5., "CDK1": 1.,
                                    "ACTB": 40. - i, "GAPDH": 10., "CD69": float(n - i), "IL2RA": 1., "FOS": 2.})
            if scores:
                # The last T cell was not scored: null must stay null.
                cell["measured_scores"] = dict(S_CYCLE=i / n, S_ACT=None if (lineage == "T_CD8" and i == n - 1) else (n - i) / n)
            result.append(cell)
    return result


def manifest():
    return dict(scores=dict(
        S_CYCLE=dict(method="test-rank-score 1.0", description="Cycling proxy", genes_scored=["MKI67", "TOP2A", "PCNA"], n_genes_set=30),
        S_ACT=dict(method="test-rank-score 1.0", description="Activation proxy", genes_scored=["CD69", "IL2RA"], n_genes_set=12,
                   caveat="Two genes only.")),
        caveat="Synthetic scores.")


def rule(action, evidence, *, lineages="*", kind="measured_score", direction="supports"):
    return dict(action=action, lineages=lineages, evidence=evidence, kind=kind, direction=direction,
                rationale="Test rule.")


def policy(*rules, **extra):
    return dict(policy_version=POLICY_VERSION, name="test_policy", rules=list(rules), **extra)


def bundle(*rules, n=12, **context_overrides):
    """Prepared bundle with measured scores and the given rules (default: two measured-score rules)."""
    rules = rules or (rule("PROLIFERATE", "S_CYCLE"), rule("ACTIVATE", "S_ACT", lineages=["T_CD8"]))
    return prepare(cells(n), context(evidence_policy=policy(*rules), **context_overrides),
                   assayed_genes=PANEL, measured_evidence=deepcopy(manifest()))


class FixedProvider:
    """In-process provider returning the same stated preference for every question."""
    name, model, max_questions = "mock_jev", "jev-1.14.0", 20
    sampling_parameters = {"temperature": None, "provider_seed": None, "status": "not exposed"}

    def __init__(self, rest=.8, prior_mode="hidden", pooling="full_distribution"):
        from virtualtissue_xenium.v2.provider import prompt_version
        self.rest, self.prior_mode, self.pooling = rest, prior_mode, pooling
        self.prompt_version = prompt_version(prior_mode, "3.0", pooling=pooling)

    def choose(self, rows, context):
        answers = {}
        for row in rows:
            events = [a for a in row["actions"] if a != "REST"]
            probabilities = {a: (self.rest if a == "REST" else (1 - self.rest) / len(events)) for a in row["actions"]}
            answers[row["cell_id"]] = dict(probabilities=probabilities, choice="REST", confidence=.5)
        return answers

    def usage(self):
        return dict(http_attempts=0, input_tokens=0, output_tokens=0, attempts_with_unknown_usage=0)
