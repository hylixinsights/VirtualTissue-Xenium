# V5 software verification — 2026-11-10

The tested v5 source includes the English sample-context questionnaire, validated optional sample details, provider forwarding and explicit lineage/action restrictions. These checks establish software behavior, not biological predictive validity.

## Local results

- **490 passed, 3 skipped** across the modern suite and a separate run of the three localhost download tests. The main run reports 487 passed, 3 skipped and 3 deselected; those three deselected tests all passed with local socket permission. The legacy preparation module was not collected because `decoupler` was unavailable in this environment.
- Historical macOS v4 equivalence: `identical: true`, zero differences for `demo270` and `object96` against the macOS-recorded snapshot. This covers prepared/control fingerprints, guidebooks, configurations, frames, expanded decisions and mocked Jev payloads for existing contexts.
- Fresh checksum-verified LUAD B5 input: 736 cells, 302 genes, 20 saved scores and all supplied polygons. Preparation reproduced 118 QC REST-only cells; expression evidence changed probabilities for 581/736 cells.
- A four-step local B5 run completed. Repeating it with questionnaire-provided descriptive sample details produced byte-identical frames at the same seed. Details were verified in the prepared/saved guidebook and constructed hidden-prior Jev payload, with the prepared fingerprint bound in the run configuration.
- Explicit restrictions remain enforced with strong cycling evidence and after beta-zero or expression-shuffle controls. Mature-neuron division remains prohibited. Mucus-related text creates no mucus action.
- Immutable engine integrity and `git diff --check` passed. Renderer, replay, ontology assets and numerical kinetic/beta tables remain unchanged from v4.

## Reproduction

After installing the documented modern extras and test dependencies:

```bash
pytest -q
python tools/v4_equivalence.py compare --same-platform-v4
python scripts/mac/check_environment.py .
```

The equivalence command and test suite execute immutable v4 commit `20098224d06fddcc00f2683194622a60febb2b25` and the candidate under the same interpreter and platform, then require exact equality of the recorded bytes and values. Full Git history is required; run `git fetch --unshallow` in a shallow clone. The fixed `tests/fixtures/v4_0_0_reference.json` is retained as a historical, environment-specific macOS artifact.

Restricted environments may need permission for the three local HTTP-server tests. See the [LUAD example](../../examples/v5/luad_tls/README.md) for real-data commands; the data are not packaged. The exact historical real-H5AD equivalence fixture was unavailable, so the fresh B5 checks above replace that check rather than claiming it passed.

CI includes v5 tests on Python 3.12/3.13, the legacy suite and macOS acceptance. Initial Linux CI on [release PR #8](https://github.com/hylixinsights/VirtualTissue-Xenium/pull/8) caught a cross-platform mismatch against the macOS snapshot. The comparison now uses the same platform for both versions and retains strict equality; no numerical logic changed. Local counts above describe the historical local run. Subsequent CI job status is recorded on that pull request.

No live Jev request was sent. The seven GO20 evidence rules remain marked `requires_biological_review`; rates, recognition and annotation mappings are declared assumptions. See [sample context](SAMPLE_CONTEXT.md) and [release notes](RELEASE_NOTES.md).
