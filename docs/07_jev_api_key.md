# The Jev API key

Cell decisions are made by **Jev**, TypeSafe's model for structured choices, through the **System One API** with the **Choice** primitive:

* API reference: <https://docs.typesafe.ai/api/system-one>
* Choice schema: <https://docs.typesafe.ai/primitives/choice>
* Pinned model versions: <https://docs.typesafe.ai/models>

Endpoint used: `POST https://api.typesafe.ai/v1/systemone`, with header `Authorization: Bearer <key>`.

## 1. Get a key

Create an account with TypeSafe ([typesafe.ai](https://typesafe.ai)) and generate an API key in their console. Check pricing and limits for your account there. VirtualTissue-Xenium only controls the **number of requests** (`MAX_HTTP_ATTEMPTS`), not money.

## 2. Store it

| Where you run | How | Name |
|---|---|---|
| Google Colab | Open your copy of **`Part2-4_Run_Jev_Simulation.ipynb`**, then 🔑 **Secrets** → *Add new secret* → enable **Notebook access** for this notebook | `TYPESAFE_API_KEY` (preferred) or `JEV_API_KEY` |
| Local Jupyter / `vtx run` | Environment variable, e.g. `export TYPESAFE_API_KEY=…` in the shell that starts Jupyter, or a `.env` loaded by your shell | same names |
| Anywhere, as a last resort | If no key is found, a **hidden prompt** (`getpass`) asks for it | — |

The lookup order is: Colab Secret `TYPESAFE_API_KEY`, Colab Secret `JEV_API_KEY`, environment `TYPESAFE_API_KEY`, environment `JEV_API_KEY`, then the prompt (`virtualtissue_xenium.colab.get_jev_key`).

Part 1, the results explorer and **Local preview (no API)** do not need a key. In Colab, granting access in Part 1 does not grant it in the simulation notebook; enable access in the notebook that will run Jev.

## 3. How the key is protected

* It is read **only inside the run cell** (or `vtx run`), passed to the provider, and cleared afterwards.
* It is **never written** to the notebook, the journal, `run_config.json`, the replay, the CSVs or the results ZIP. If an API response echoes the key, it is replaced by `[redacted]` before anything is saved. Error messages contain only status codes. The tests scan every output for a fake key.
* Requests follow **no redirects**, so the `Authorization` header cannot be forwarded to another host.
* `.gitignore` excludes `.env`, `*.key` and `secrets*.json`. Never commit a key, and never paste one into a code cell or form field: notebook outputs are saved in Drive and are often shared.
* If a key leaks, revoke it in the TypeSafe console and create a new one.

## 4. What is sent

Each request carries up to 20 independent cell questions (≤ 40,000 bytes):

* `model`: the pinned version you chose (`jev-1.14.0` by default);
* `state`: two fixed sentences describing the setting and what the evidence means;
* per cell: the permitted options with one-line descriptions, plus its local observation. The observation contains the type label, provisional identity, 20 t0 AUCell percentiles and raw values, simulated state, simulated fields, contacts (IDs and labels of neighbours), and which options were gated out and why.

**Not sent:** raw counts, gene lists, coordinates, images, sample metadata, the contract or rulebook, or any desired outcome. Cell IDs are internal integers; the mapping to original Xenium IDs stays local. Before you run on unpublished data, review `first_request_preview.json` from the preflight. It is exactly what the first request will contain.

## 5. Model versions

Pin a version (`jev-X.Y.Z`). The run stops (`provider_failed`) if the API answers with another model, so a run never mixes model versions. `jev-1.14.0` is the default; `jev-1.13.0` was used in the original LUAD3 prototype. A different model is a different configuration: use a new `RUN_NAME`.

## 6. Usage reporting

`run_report.json → usage` reports `http_attempts` (requests sent, including failures), `reused_batches` (answers taken from the journal on resume, which cost nothing), and `input_tokens` / `output_tokens` summed from the API's `usage` fields when present.
