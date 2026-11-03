# Synthetic object with polygons

This fixture contains **96 generated cells, not experimental data**. Counts,
labels and polygons are arbitrary engineering inputs. They do not validate
biological activity, cell identity, tissue structure or future trajectories.

Run from the repository root, with the project installed or `PYTHONPATH=src`:

```sh
PYTHONPATH=src python examples/v3/make_object_demo.py \
  --output outputs/object-demo/object.json \
  --context-output outputs/object-demo/context.json \
  --include-invalid-polygon

PYTHONPATH=src python -m virtualtissue_xenium.v2.cli roi-catalog \
  outputs/object-demo/object.json outputs/object-demo/catalog.html \
  --cells-per-roi 32 --grid-um 40 --max-candidates 12
```

The generator refuses to overwrite outputs. The default seed is `20261101`;
`--seed` changes the generated counts and geometry. All coordinates and polygon
vertices use micrometres. The optional context file is an explicit, synthetic
adaptation of `context_user_object.json`; the same context is also recorded in
the object JSON. Source metadata includes the generator version, seed and hash
of the packaged gene signatures.

Composition is fixed: 40 malignant, 18 T_CD8, 12 myeloid, 8 fibroblast,
6 endothelial, 4 T_reg, 4 NK and 4 B. The three artificial patches have different
compositions and do not encode biological regions. Existing v3 lineage names
are assigned before expression is generated; they are not inferred from it.

The panel includes the packaged human signature genes plus explicitly named
`SYNTHETIC_BACKGROUND_*` features. Heterogeneous counts exercise whole-type
empirical references for the abundant type and the weak within-cell panel
heuristic for rare types. One cell deliberately has only three counts and fails
the context's 20-count threshold. Reference support changes after selecting a
smaller region; the catalog does not fit scores for the whole raw object.

All rings are closed. With `--include-invalid-polygon`, one ring is deliberately
self-intersecting, identified in cell and annotation metadata, and should use
centroid fallback. Without that option, all 96 rings are valid synthetic shapes.

Open the generated catalog and inspect compositions before choosing a region.
**No ROI is selected automatically.** After saving a selection in the catalog,
pass that selection file and the explicit context to `select-roi`; its output
can then be inspected and replayed using the normal v3 commands. The generator
does not call Jev or run a simulation.
