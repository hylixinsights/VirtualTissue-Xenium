# VirtualTissue-Xenium 4.0.0

Released 2026-11-06. V4 is the current primary software release. Use `vtx4` for the modern workflow; `vtx2` remains a compatible alias and `vtx` retains the legacy physical engine.

The replay restores the original rounded cell bodies, neighboring-cell clipping, visible nuclei, tissue-scale sizes, context transparency and event glows, together with the dark dashboard and full sidebar. **Original cells** transitions supplied t0 segmentation into illustrative simulated silhouettes. **Provided outlines** shows exact saved polygons; **Centroids** shows recorded positions. These presentation modes preserve the recorded coordinates, supplied geometry and interaction graph.

The v3 simulation mechanisms and corrections remain in place. Simulation schema **3.0**, ontologies, kinetic assumptions and existing `docs/v3` and `examples/v3` paths are intentionally unchanged. V4 identifies the software distribution; it does not establish biological predictive validity. Kinetic defaults remain **DRAFT assumptions**. See the [rate review](../v3/RATE_REVIEW.md) and [limitations](../../LIMITATIONS.md).

Refresh a saved v3 trajectory after updating the package:

```bash
vtx4 replay outputs/v3/ensemble_990_1350_48h/runs/recognition/seed_4
```

This reads saved frames and the report and regenerates presentation files. It makes no Jev requests and changes no decisions, readouts or simulation exports. Keep the generated HTML with its `replay_frames` folder. The browser retains the bounded two-chunk frame cache.

Jev event timing and the volume or duplication of scientific exports are unchanged in this release.
