# Jev decisions: what the model is told, how its answer is used, how to measure it

## 1. The request

A bundle without an evidence policy or optional sample-context fields sends the v4 request (`vtx-context-guidebook-3.1-evidence-basis`), byte for byte. Declared `sample_context` fields are forwarded with their investigator-declaration semantics in either evidence format; see [sample context](SAMPLE_CONTEXT.md). A bundle with a policy sends prompt `vtx-context-guidebook-5.0-evidence-support`:

| Where | Content |
|---|---|
| `state.evidence_definitions` | Once per request: for each item, `kind` (`gene_set` or `measured_score`), the assayed genes that contributed, the size of the source set, and for measured scores the description and any caveat. |
| `state.measured_scores` | The scoring method and the manifest's general caveat. |
| `state.evidence_unavailable` | Items no question in the request can use, with the reason, stated once. |
| `observation.t0_evidence` | Per cell: `name → z`, two decimals, available items only. |
| `observation.evidence_support` | Per cell: the reference `basis` and `n_ref` shared by its values, `n_detected` for gene sets, and `other` for items with a different reference. |
| `observation.evidence_unavailable` | Per cell: items missing for this cell only, with the reason. |
| `state.evidence_semantics` | States that t0 evidence was measured once and is never updated, that `state` and `local` are produced by the simulator, that a small `n_ref` or `n_detected` means weak evidence, and that unavailable is unknown, not zero. |

The evidence policy's rules are not sent. `provider_evidence` in the policy limits which items are sent; fewer items cost fewer tokens.

`vtx5 cost prepared.json --steps N` estimates calls and input tokens from the t0 questions (add `--jev-pooling event_share` to count that rule's extra instruction). On a 736-cell region with 23 available evidence items the estimate is about 1,900 input tokens per question.

## 2. How the answer is used

Jev returns preferences over the offered menu. The run combines them with the local prior.

**`full_distribution` (default, the v4 rule).** `P(a) ∝ prior(a)^(1−α) × jev(a)^α` over every action including REST. With the prior hidden, α is 1 and Jev's answer is used as given. Jev's REST preference therefore sets how often events happen. The local prior scales with the time step; a stated preference does not, so the event rate changes when the time step changes. With a fixed answer of REST 0.8 / ACTIVATE 0.2 and a local hazard of 0.1 per hour:

| α | dt 0.25 h | dt 0.5 h | dt 1 h |
|---|---:|---:|---:|
| 0 (local) | 0.099 | 0.098 | 0.095 |
| 0.5 | 0.295 | 0.203 | 0.140 |
| 1 (hidden prior) | 0.800 | 0.400 | 0.200 |

Values are proposed events per hour.

**`event_share` (opt-in).** The kernel keeps its event probability, `P(event) = 1 − exp(−Σk·dt)`, and Jev's preference is pooled only among the non-REST actions:

```text
P(REST) = prior(REST)
P(a)    = (1 − prior(REST)) × share(a)
share(a) ∝ prior_share(a)^(1−α) × jev_share(a)^α
```

Jev then decides which event occurs, not whether one occurs. In the table above every row becomes the local row. If Jev puts no preference on any event, the prior shares are kept. With the prior hidden, α is 1 and the composition of events is entirely Jev's; the timing still comes from the kinetic table, which is itself an assumption.

```bash
vtx5 run prepared.json out --mode jev --jev-pooling event_share --max-http-attempts 500
```

In an ensemble manifest, add `"jev_pooling": "event_share"` to a Jev arm. Compare the two rules as arms of one ensemble; `compare` treats the pooling rule as a run setting and refuses to mix runs made with different rules. The prompt tells the model which rule applies, the prompt version gains `-event-share`, and `config.json` records the rule; a default run records nothing new and keeps its v4 configuration hash.

**Every Jev run reports the difference.** `report.json → jev_event_mass.by_lineage` gives the mean probability of doing anything but REST as declared by the kernel (`prior`), stated by Jev (`jev`) and applied, with `applied_over_prior` and both as events per hour. Under `event_share` the ratio is 1 by construction.

## 3. Probing the model on a frozen state

A trajectory cannot attribute Jev's behaviour: after one step two arms are in different states with different menus. `jev-probe` asks about one frozen t0 state instead. The same cells, menus and request batches are used in every arm.

| Arm | What changes |
|---|---|
| `baseline` | Nothing. |
| `evidence_masked` | Every t0 evidence value is withheld. With an evidence policy the definitions go with them, so the arm compares evidence present against evidence absent. |
| `evidence_exchanged` | The cell receives the evidence of another cell of the same declared type. |
| `neighbourhood_exchanged` | The cell receives the local summary of another cell of the same lineage and menu. |
| `context_exchanged` | Tissue, disease and perturbation text come from `--donor-context`. |

```bash
vtx5 jev-probe plan probe_dir --prepared prepared.json --cells-per-lineage 6 --repeats 3
# plan.json lists cells, batches, request hashes and the estimate; payloads.json holds every request
export JEV_API_KEY=...
vtx5 jev-probe run probe_dir --max-http-attempts 120
vtx5 jev-probe analyze probe_dir
```

- `plan` and `analyze` send nothing. `run` refuses a budget below the whole plan and never sends a request twice, so an interrupted run can be continued.
- A request that fails is not resent, because its usage is unknown. `run` reports it and exits with status 2 after sending the rest. `analyze` then leaves that request group out of every arm and lists it under `request_groups_excluded`.
- An exchange arm skips a cell that has no counterpart, or whose counterpart shows the same content; `plan.json` lists these under `cells_unchanged_by_design`. On a panel where no evidence is usable, `evidence_exchanged` changes nothing and reports zero cells.
- Repeats go to independent journals. Use at least two; with one, replicate noise is unmeasured.
- With the prior shown, the guidebook prior stays the baseline prior in every arm, so an evidence arm isolates the evidence shown directly. Use the hidden prior to see the model's own preferences.

`analysis.json` reports, for each prior mode:

- **`replicate_noise`**: mean Jensen-Shannon divergence between identical baseline requests, and the share of repeat pairs with the same most-preferred action.
- **Per arm**: mean JS between the arm and the baseline; `js_excess_over_replicates`, the per-cell difference between that and the average divergence between repeats within the baseline and within the arm, with a bootstrap interval that resamples request groups, since questions sent together are not independent; the share of arm-baseline repeat pairs with the same most-preferred action; mean event probability, with the baseline's over the same cells.
- **`evidence_tracking`**: per lineage and action with declared evidence, the rank correlation across probed cells between the evidence z, signed by the rule's direction, and the share of event preference Jev gives the action. The same quantity for the local prior is listed beside it. `evidence_shown_to_provider` is false when the policy's `provider_evidence` withheld that item from Jev.

Read an arm's excess against zero. An excess near zero means the manipulation did not move the answers beyond replicate noise on this question bank. The probe states no threshold and certifies nothing: it describes one model version and one set of cells at t0.

`jev-audit` (v4) still summarises repeated-payload variance from run journals.

## 4. What was and was not tested

The prompt, both pooling rules and the probe were tested with in-process providers that return fixed or evidence-following answers. No live Jev request was made while preparing 5.0.0. How the real model responds to the new prompt, and whether it responds to evidence at all, is what the probe is for.
