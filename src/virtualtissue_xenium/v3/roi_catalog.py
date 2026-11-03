"""Offline ROI composition catalogues and explicitly user-declared selection.

Labels and polygons are supplied data, not predictions or certified annotations.
Candidates are geometric alternatives; none is selected or called an invasive
front automatically. Scores shown here describe the source prepared bundle.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from copy import deepcopy
import json
import heapq
import math
from pathlib import Path
from statistics import median

from ..v2.inputs import digest, read_json, validate_cells, validate_prepared
from ..v2.scoring import depth_bins

CATALOG_VERSION = "vtx-roi-catalog-1"
SELECTION_VERSION = "vtx-roi-selection-1"


def _source(bundle):
    bundle = read_json(bundle) if isinstance(bundle, (str, Path)) else bundle
    if not isinstance(bundle, dict):
        raise ValueError("Use a raw JSON cells object or a schema 3.0 prepared bundle.")
    if "guidebook" in bundle or "fingerprint" in bundle:
        validate_prepared(bundle)
        if bundle["schema_version"] != "3.0":
            raise ValueError("Prepared ROI catalogues require schema 3.0.")
        if bundle["context"]["interaction_mode"] != "spatial":
            raise ValueError("ROI selection requires measured spatial coordinates.")
        bundle = dict(bundle, source_kind="prepared_bundle", score_coverage_status="computed_in_source_bundle",
                      qc_interpretation="source_bundle_qc")
    else:
        # Inspect counts and shape without fitting references or constructing a
        # guidebook over the entire section. The raw input remains hash-bound.
        from .geometry import normalize_geometry
        cells, _ = normalize_geometry(validate_cells(bundle.get("cells")))
        diagnostic = dict(min_transcripts=20, depth_bins=10, min_reference_cells=30,
                          interpretation="diagnostic defaults only; actual scoring is determined when the selected ROI is prepared")
        entries = {c["cell_id"]: dict(expression_total=sum(c["expression"].values()),
                     detected_genes=sum(v > 0 for v in c["expression"].values()),
                     qc="pass" if sum(c["expression"].values()) >= diagnostic["min_transcripts"] else "below_min_transcripts")
                   for c in cells}
        if any(not math.isfinite(entry["expression_total"]) for entry in entries.values()):
            raise ValueError("Transcript totals overflowed; inspect the raw counts.")
        # Invalid polygon numerics are preserved as inert raw_input_json by the
        # geometry adapter. Hash that normalized representation so its declared
        # fallback remains available without admitting nonfinite JSON metadata.
        try:
            source_hash = digest(dict(bundle, cells=cells))
        except (ValueError, TypeError) as exc:
            raise ValueError("Raw metadata outside normalized geometry must be finite JSON values.") from exc
        bundle = dict(cells=cells, fingerprint=source_hash, source_kind="raw_json",
            context={"scoring": diagnostic}, guidebook={"cells": entries, "provenance": {"gene_sets": {}}},
            measurement=deepcopy(bundle.get("measurement", {"unit": "not_declared"})),
            annotation_provenance=deepcopy(bundle.get("annotation_provenance", {"origin": "user_supplied"})),
            score_coverage_status="not_computed", qc_interpretation="diagnostic_transcript_threshold_not_prepared_qc")
    from .geometry import validate_geometry_records
    validate_geometry_records(bundle["cells"])
    for cell in bundle["cells"]:
        if any(type(cell.get(k)) not in (int, float) or not math.isfinite(cell[k]) for k in ("x", "y")):
            raise ValueError("Every ROI catalogue cell needs finite measured x/y in micrometres.")
        _polygon(cell)
    return bundle


def _polygon(cell):
    if "t0_geometry" in cell:
        from .geometry import projected_polygon
        return projected_polygon(cell)
    polygon = cell.get("polygon")
    if polygon is None:
        return None
    if (not isinstance(polygon, list) or len(polygon) < 3 or
            any(not isinstance(p, (list, tuple)) or len(p) != 2 or
                any(type(v) not in (int, float) or not math.isfinite(v) for v in p) for p in polygon)):
        raise ValueError("polygon must be at least three finite [x,y] vertices in micrometres.")
    points = [[float(x), float(y)] for x, y in polygon]
    if len({tuple(p) for p in points}) < 3:
        raise ValueError("polygon needs at least three distinct vertices.")
    return points


def _geometry_kind(cell):
    if "t0_geometry" in cell:
        return "supplied_polygon" if cell["t0_geometry"]["status"] == "valid" else "invalid_polygon_centroid_fallback"
    return "supplied_polygon" if _polygon(cell) is not None else "centroid_fallback"


def _parameters(cells_per_roi, grid_um, max_candidates):
    if type(cells_per_roi) is not int or cells_per_roi < 1:
        raise ValueError("cells_per_roi must be a positive integer.")
    if type(max_candidates) is not int or not 1 <= max_candidates <= 1000:
        raise ValueError("max_candidates must be an integer from 1 to 1000 per sample.")
    if type(grid_um) not in (int, float) or not math.isfinite(grid_um) or grid_um <= 0:
        raise ValueError("grid_um must be positive and finite.")
    return dict(cells_per_roi=cells_per_roi, grid_um=float(grid_um), max_candidates=max_candidates)


def _geometry_bounds(cells):
    points = [(c["x"], c["y"]) for c in cells]
    for c in cells:
        points.extend(_polygon(c) or [])
    return [min(p[0] for p in points), min(p[1] for p in points),
            max(p[0] for p in points), max(p[1] for p in points)]


def _candidates(cells, parameters):
    by_sample = defaultdict(list)
    for cell in cells:
        by_sample[cell["sample_id"]].append(cell)
    grid = parameters["grid_um"]
    for sample, members in sorted(by_sample.items()):
        centers = sorted({(math.floor(c["x"] / grid + .5) * grid,
                           math.floor(c["y"] / grid + .5) * grid) for c in members})
        cap = parameters["max_candidates"]
        if len(centers) > cap:
            indices = [0] if cap == 1 else [i * (len(centers) - 1) // (cap - 1) for i in range(cap)]
            centers = [centers[i] for i in indices]
        for center in centers:
            chosen = heapq.nsmallest(parameters["cells_per_roi"], members,
                key=lambda c: ((c["x"]-center[0])**2 + (c["y"]-center[1])**2, c["cell_id"]))
            candidate_id = "roi_" + digest([sample, center, parameters])[:20]
            yield candidate_id, sample, list(center), sorted(chosen, key=lambda c: c["cell_id"])


def _summary(cells, bundle):
    entries = bundle["guidebook"]["cells"]
    scoring = bundle["context"]["scoring"]
    n = len(cells)
    by_type = defaultdict(list)
    for cell in cells:
        by_type[cell["cell_type"]].append(cell)
    eligible = [c for c in cells if entries[c["cell_id"]]["qc"] == "pass"]
    bins = depth_bins(eligible, scoring["depth_bins"])
    potential = []
    for name, members in sorted(by_type.items()):
        counts = Counter(bins[c["cell_id"]] for c in members if c["cell_id"] in bins)
        potential.append(dict(cell_type=name, n_units=len(members), n_qc_pass=sum(counts.values()),
            within_roi_depth_bin_counts={str(k): v for k, v in sorted(counts.items())},
            cells_in_strata_meeting_minimum=sum(v for v in counts.values() if v >= scoring["min_reference_cells"]),
            minimum_reference_cells=scoring["min_reference_cells"]))
    signatures = bundle["guidebook"]["provenance"]["gene_sets"]
    coverage = {}
    for signature in sorted(signatures):
        evidence = [entries[c["cell_id"]]["evidence"][signature] for c in cells]
        coverage[signature] = dict(available=sum(e.get("status") in ("available", "available_panel_heuristic") for e in evidence),
            empirical_available=sum(e.get("status") == "available" for e in evidence),
            weak_panel_available=sum(e.get("status") == "available_panel_heuristic" for e in evidence),
            informative=sum(bool(e.get("informative")) for e in evidence),
            statuses=dict(sorted(Counter(e.get("status", "unknown") for e in evidence).items())),
            source_reference_size_min=min(e.get("reference_size", 0) for e in evidence),
            source_reference_size_max=max(e.get("reference_size", 0) for e in evidence),
            panel_coverage=min(e.get("panel_coverage", 0) for e in evidence))
    totals = [float(entries[c["cell_id"]]["expression_total"]) for c in cells]
    detected = [entries[c["cell_id"]]["detected_genes"] for c in cells]
    unknown = sum(c["lineage"] == "unknown" for c in cells)
    lineage_counts = Counter(c["lineage"] for c in cells)
    return dict(n_units=n, by_lineage={key: dict(count=count, fraction=count/n) for key, count in sorted(lineage_counts.items())},
        unknown_count=unknown, unknown_fraction=unknown/n, qc_pass=len(eligible), qc_fail=n-len(eligible),
        transcript_total=sum(totals), transcript_median=median(totals), detected_genes_median=median(detected),
        polygons=sum(_geometry_kind(c) == "supplied_polygon" for c in cells),
        centroid_fallback=sum(_geometry_kind(c) != "supplied_polygon" for c in cells),
        invalid_polygons=sum(_geometry_kind(c) == "invalid_polygon_centroid_fallback" for c in cells),
        source_score_coverage=coverage, score_coverage_status=bundle["score_coverage_status"],
        qc_interpretation=bundle["qc_interpretation"], min_transcripts=scoring["min_transcripts"],
        potential_references=potential,
        reference_note="Potential counts are descriptive, not an availability gate for every reference strategy. Scores shown are frozen source-bundle scores; subsetting and preparing again may change reference bins, scores and coverage.")


def _descriptor(bundle, candidate_id, sample, center, chosen, parameters, data_hash, summary=None):
    return dict(schema_version=SELECTION_VERSION, source_fingerprint=bundle["fingerprint"],
        source_kind=bundle["source_kind"], data_sha256=data_hash, candidate_id=candidate_id, sample_id=sample,
        selection_method="user_selected_geometric_candidate", selection_declared_by="user",
        generator=dict(name="occupied_grid_nearest_units", version="1", parameters=parameters),
        center_um=center, cell_ids=[c["cell_id"] for c in chosen],
        window=dict(kind="selected_geometry_bounds", bounds_um=_geometry_bounds(chosen),
                    is_statistical_observation_window=False),
        diagnostics=deepcopy(summary if summary is not None else _summary(chosen, bundle)),
        geometry_rule="Membership uses measured centroids, with cell_id distance ties; polygons are displayed without changing membership.")


def build_roi_catalog(bundle, *, cells_per_roi=1000, grid_um=500, max_candidates=50):
    """Describe deterministic alternatives; retain all supplied labels, even unknown.

If a sample has fewer than cells_per_roi, a candidate includes its available
units and is marked undersized. No cell identity is inferred and no candidate
is chosen for the user. The complete source is never changed.
"""
    bundle = _source(bundle)
    data_hash = digest(bundle["cells"])
    params = _parameters(cells_per_roi, grid_um, max_candidates)
    samples, map_cells = {}, []
    by_sample = defaultdict(list)
    for cell in sorted(bundle["cells"], key=lambda c: c["cell_id"]):
        by_sample[cell["sample_id"]].append(cell)
        entry = bundle["guidebook"]["cells"][cell["cell_id"]]
        map_cells.append(dict(cell_id=cell["cell_id"], sample_id=cell["sample_id"],
            cell_type=cell["cell_type"], lineage=cell["lineage"], x=cell["x"], y=cell["y"],
            polygon=_polygon(cell), geometry=_geometry_kind(cell),
            geometry_warnings=deepcopy(cell.get("t0_geometry", {}).get("warnings", [])),
            area_um2=cell.get("t0_geometry", {}).get("area_um2"),
            qc=entry["qc"], transcript_total=entry["expression_total"], detected_genes=entry["detected_genes"]))
    for sample, members in sorted(by_sample.items()):
        samples[sample] = dict(bounds_um=_geometry_bounds(members), summary=_summary(members, bundle))
    candidates = []
    for cid, sample, center, chosen in _candidates(bundle["cells"], params):
        summary = _summary(chosen, bundle)
        candidates.append(dict(candidate_id=cid, sample_id=sample, center_um=center,
            undersized=len(chosen) < cells_per_roi, summary=summary,
            selection=_descriptor(bundle, cid, sample, center, chosen, params, data_hash, summary)))
    result = dict(schema_version=CATALOG_VERSION, source_fingerprint=bundle["fingerprint"],
        source_kind=bundle["source_kind"], data_sha256=data_hash, parameters=params,
        measurement=deepcopy(bundle["measurement"]), score_coverage_status=bundle["score_coverage_status"],
        annotation_provenance=deepcopy(bundle.get("annotation_provenance", {})),
        annotation_interpretation="User-supplied labels and geometry; no inferred malignancy or verification of annotation provenance.",
        scoring=deepcopy(bundle["context"]["scoring"]), selected_candidate=None,
        samples=samples, cells=map_cells, candidates=candidates)
    result["catalog_sha256"] = digest(result)
    return result


def validate_roi_selection(bundle, selection):
    """Validate source binding, exact IDs and deterministic candidate geometry.

This verifies the declaration, not whether a person actually reviewed a file.
The returned descriptor is copied. It never modifies or prepares the source.
"""
    bundle = _source(bundle)
    return _validate_selection(bundle, selection)


def _validate_selection(bundle, selection):
    selection = read_json(selection) if isinstance(selection, (str, Path)) else selection
    if not isinstance(selection, dict) or selection.get("schema_version") != SELECTION_VERSION:
        raise ValueError("Expected a vtx-roi-selection-1 selection descriptor.")
    data_hash = digest(bundle["cells"])
    if selection.get("source_fingerprint") != bundle["fingerprint"] or selection.get("data_sha256") != data_hash:
        raise ValueError("ROI selection refers to a different source fingerprint or data hash.")
    generator = selection.get("generator")
    if not isinstance(generator, dict) or generator.get("name") != "occupied_grid_nearest_units" or generator.get("version") != "1":
        raise ValueError("Unrecognized ROI candidate generator.")
    params = generator.get("parameters")
    if not isinstance(params, dict) or set(params) != {"cells_per_roi", "grid_um", "max_candidates"}:
        raise ValueError("Selection needs the original catalogue parameters.")
    params = _parameters(**params)
    ids = selection.get("cell_ids")
    if not isinstance(ids, list) or not ids or any(not isinstance(cid, str) for cid in ids) or len(set(ids)) != len(ids):
        raise ValueError("Selection cell_ids must be nonempty, unique exact source IDs.")
    for cid, sample, center, chosen in _candidates(bundle["cells"], params):
        if cid != selection.get("candidate_id"):
            continue
        expected = _descriptor(bundle, cid, sample, center, chosen, params, data_hash)
        if selection != expected:
            raise ValueError("Selection IDs, window, sample or declaration differ from the recorded geometric candidate.")
        return deepcopy(expected)
    raise ValueError("Selection candidate is absent from the source catalogue.")


def extract_roi_selection(bundle, selection):
    """Return independent copies of selected cells, ordered by exact source ID."""
    bundle = _source(bundle)
    descriptor = _validate_selection(bundle, selection)
    by_id = {c["cell_id"]: c for c in bundle["cells"]}
    return [deepcopy(by_id[cid]) for cid in descriptor["cell_ids"]]


def write_roi_catalog(bundle, output_html, **options):
    """Write a standalone HTML catalogue and its JSON; no service calls."""
    catalog = build_roi_catalog(bundle, **options)
    output = Path(output_html)
    if output.suffix.lower() != ".html":
        raise ValueError("Choose an .html output path for the offline ROI catalogue.")
    sidecar = output.with_suffix(".catalog.json")
    if output.exists() or sidecar.exists():
        raise ValueError("ROI catalogue output already exists; choose new paths.")
    output.parent.mkdir(parents=True, exist_ok=True)
    encoded = json.dumps(catalog, ensure_ascii=False, allow_nan=False).replace("<", "\\u003c")
    with sidecar.open("x", encoding="utf8") as stream:
        json.dump(catalog, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")
    with output.open("x", encoding="utf8") as stream:
        stream.write(_HTML.replace("__CATALOG_JSON__", encoded))
    return dict(html=str(output.resolve()), catalog=str(sidecar.resolve()),
                catalog_sha256=catalog["catalog_sha256"], candidates=len(catalog["candidates"]),
                selected_candidate=None)


_HTML = r'''<!doctype html>
<html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>VirtualTissue · Region selection</title>
<style>
*{box-sizing:border-box}body{margin:0;font:15px system-ui,sans-serif;background:#edf3f1;color:#16332f}header{padding:22px 28px;background:#123f3a;color:white}h1{font-size:24px;margin:0 0 8px}p{line-height:1.5;margin:8px 0}main{display:grid;grid-template-columns:minmax(380px,1.6fr) minmax(320px,1fr);gap:18px;padding:20px;max-width:1600px;margin:auto}.panel{background:white;padding:18px;border:1px solid #cedbd6;border-radius:12px;min-width:0}h2{font-size:18px;margin:0 0 12px}button,select{font:inherit;padding:9px 12px;border:1px solid #91aaa2;border-radius:7px;background:white;color:#123f3a}button{cursor:pointer}button.primary{background:#147b69;color:white;border-color:#147b69}button:disabled{opacity:.45;cursor:default}label{display:inline-flex;align-items:center;gap:5px;margin:4px 8px 4px 0}canvas{display:block;width:100%;height:530px;background:#f8faf9;border:1px solid #d4dfda;margin:12px 0;border-radius:8px}small,.muted{color:#526b62;font-size:12px}.legend{display:flex;gap:10px;flex-wrap:wrap;margin:8px 0}.dot{display:inline-block;width:9px;height:9px;border-radius:50%;margin-right:5px}.candidate{display:block;width:100%;text-align:left;margin:6px 0;background:#f6faf8}.candidate.active{border:2px solid #147b69;background:#e4f5ee}.candidate span{display:block;font-size:12px;margin-top:4px}#choices{max-height:420px;overflow:auto}.metrics{display:flex;flex-wrap:wrap;gap:8px}.metric{background:#eef5f2;border-radius:6px;padding:8px}.metric b{display:block;font-size:20px}.notice{padding:9px 12px;background:#fff4d9;border-radius:6px;font-size:13px}table{border-collapse:collapse;width:100%;font-size:12px;margin:12px 0}td,th{border-bottom:1px solid #e1e9e5;padding:7px 4px;text-align:left}th{font-weight:600}.scroll{overflow:auto}details{margin:12px 0}pre{white-space:pre-wrap;overflow-wrap:anywhere;font-size:11px}#details{grid-column:1/-1}#pointer{min-height:20px} @media(max-width:900px){main{grid-template-columns:1fr;padding:10px}canvas{height:420px}#details{grid-column:auto}}
</style>
<header><h1>Choose a region to simulate</h1><p>Compare composition and coverage before selecting. No region has been selected automatically.</p></header>
<main><section class="panel"><h2>Sample map</h2><label>Sample <select id="sample"></select></label><label><input id="polygons" type="checkbox" checked> Supplied polygons</label><label><input id="zoom" type="checkbox"> Zoom to region</label><button id="clear">Clear selection</button><div id="overview" class="muted"></div><canvas id="map" aria-label="Map of cells and candidate centers"></canvas><div id="pointer" class="muted"></div><div id="legend" class="legend"></div><p class="muted">Numbered centers are geometric alternatives. Click a center or an option alongside the map. When no polygon is available, the point represents the measured centroid.</p></section>
<section class="panel"><h2>Candidate regions</h2><p class="muted">Fixed spatial order; no biological ranking.</p><div id="choices"></div><p id="selection-status">Choose a region to enable export.</p><button class="primary" id="download" disabled>Download selection JSON</button><p class="muted">The file contains exact IDs and source hashes. Selecting a region does not modify your data or run a simulation.</p></section>
<section class="panel" id="details"><h2 id="detail-title">Sample summary</h2><p id="evidence-status" class="notice"></p><div id="metrics" class="metrics"></div><div id="composition" class="scroll"></div><p class="notice">Existing object labels and supplied polygons are preserved. The catalog does not infer malignancy, certify annotations or require H&amp;E or independent review. The displayed rectangle is a visual boundary, not a validated statistical window.</p><details open><summary>Signature coverage in the source</summary><p class="muted">Scores from the source prepared bundle. Preparing only the region may change references, coverage and scores. “Informative” means there is variation in the reference, not biological validation. “Weak heuristic” uses the gene panel and does not represent an empirical reference with more cells.</p><div id="coverage" class="scroll"></div></details><details><summary>Potential references within the region</summary><p class="muted">Counts after QC by type and depth bin. These describe reference size; alternative strategies for rare types may use other references.</p><div id="references" class="scroll"></div></details><details><summary>Provenance and selection</summary><pre id="provenance"></pre></details></section></main>
<script id="catalog-data" type="application/json">__CATALOG_JSON__</script>
<script>
'use strict';
const catalog=JSON.parse(document.getElementById('catalog-data').textContent), el=id=>document.getElementById(id);
const canvas=el('map'), ctx=canvas.getContext('2d'), sample=el('sample'), download=el('download');
let selected=null, candidatePoints=[], cellPoints=[], view=null;
const lineages=[...new Set(catalog.cells.map(c=>c.lineage))].sort(), colors=Object.fromEntries(lineages.map((x,i)=>[x,x==='unknown'?'#89938e':`hsl(${(i*137.508)%360},52%,42%)`]));
for(const name of Object.keys(catalog.samples)){const option=document.createElement('option');option.value=name;option.textContent=name;sample.append(option)}
for(const name of lineages){const item=document.createElement('span'),dot=document.createElement('i');dot.className='dot';dot.style.background=colors[name];item.append(dot,document.createTextNode(name));el('legend').append(item)}
function percent(x){return (100*x).toFixed(1)+'%'}
function candidates(){return catalog.candidates.filter(c=>c.sample_id===sample.value)}
function units(){return catalog.cells.filter(c=>c.sample_id===sample.value)}
function table(target,headers,rows){const host=el(target);host.replaceChildren();const table=document.createElement('table'),head=document.createElement('tr');for(const h of headers){const th=document.createElement('th');th.textContent=h;head.append(th)}table.append(head);for(const row of rows){const tr=document.createElement('tr');for(const value of row){const td=document.createElement('td');td.textContent=String(value);tr.append(td)}table.append(tr)}host.append(table)}
function details(){const summary=selected?selected.summary:catalog.samples[sample.value].summary;el('evidence-status').textContent=summary.score_coverage_status==='not_computed'?`Coverage has not been calculated yet. This catalog uses only counts and geometry; the ${summary.min_transcripts}-transcript threshold is descriptive screening. Scores will be calculated when the chosen region is prepared.`:'Coverage was calculated in the prepared source; selection does not recalculate or modify scores.';el('detail-title').textContent=selected?'Selected region · '+selected.candidate_id:'Sample summary';el('metrics').replaceChildren();for(const [label,value] of [['Units',summary.n_units],['Unknown',percent(summary.unknown_fraction)],[summary.score_coverage_status==='not_computed'?'Above descriptive threshold':'Pass QC',summary.qc_pass],['Transcripts (median)',summary.transcript_median],['Detected genes (median)',summary.detected_genes_median],['Polygons / centroids',summary.polygons+' / '+summary.centroid_fallback],['Invalid polygons',summary.invalid_polygons]]){const item=document.createElement('div'),b=document.createElement('b');item.className='metric';b.textContent=value;item.append(b,document.createTextNode(label));el('metrics').append(item)}table('composition',['Lineage','Units','Fraction'],Object.entries(summary.by_lineage).map(([k,v])=>[k,v.count,percent(v.fraction)]));table('coverage',['Signature','Available','Empirical','Weak heuristic','Informative','Panel','Source reference'],Object.entries(summary.source_score_coverage).map(([k,v])=>[k,`${v.available}/${summary.n_units}`,v.empirical_available,v.weak_panel_available,`${v.informative}/${summary.n_units}`,percent(v.panel_coverage),`${v.source_reference_size_min}–${v.source_reference_size_max}`]));table('references',['Type','After QC','Counts per bin','Units in bins ≥ minimum'],summary.potential_references.map(r=>[r.cell_type,r.n_qc_pass,Object.entries(r.within_roi_depth_bin_counts).map(([k,v])=>`${k}: ${v}`).join('; ')||'—',r.cells_in_strata_meeting_minimum+' (min. '+r.minimum_reference_cells+')']));el('provenance').textContent=JSON.stringify({source_kind:catalog.source_kind,source_fingerprint:catalog.source_fingerprint,data_sha256:catalog.data_sha256,annotation_provenance:catalog.annotation_provenance,measurement:catalog.measurement,selection:selected?selected.selection:null},null,2);download.disabled=!selected;el('selection-status').textContent=selected?`${summary.n_units} units selected. ${selected.undersized?'The sample has fewer units than the requested size.':''}`:'Choose a region to enable export.'}
function choose(candidate){selected=candidate;for(const button of el('choices').children)button.classList.toggle('active',selected!==null&&button.dataset.candidate===selected.candidate_id);details();draw()}
function list(){el('choices').replaceChildren();candidates().forEach((candidate,i)=>{const button=document.createElement('button'),info=document.createElement('span'),s=candidate.summary;button.className='candidate';button.dataset.candidate=candidate.candidate_id;button.append(document.createTextNode(`${i+1}. Center (${candidate.center_um.map(v=>v.toFixed(0)).join(', ')}) µm · ${s.n_units} units`));info.textContent=`Unknown ${percent(s.unknown_fraction)} · QC ${s.qc_pass}/${s.n_units} · `+Object.entries(s.by_lineage).map(([k,v])=>`${k} ${v.count}`).join(', ');button.append(info);button.onclick=()=>choose(candidate);el('choices').append(button)});const s=catalog.samples[sample.value].summary;el('overview').textContent=`${s.n_units} units · ${s.polygons} supplied polygons · ${s.centroid_fallback} centroids without polygons`;choose(null)}
function draw(){const box=canvas.getBoundingClientRect(),dpr=window.devicePixelRatio||1;canvas.width=box.width*dpr;canvas.height=box.height*dpr;ctx.setTransform(dpr,0,0,dpr,0,0);ctx.clearRect(0,0,box.width,box.height);const raw=selected&&el('zoom').checked?selected.selection.window.bounds_um:catalog.samples[sample.value].bounds_um;const [xmin,ymin,xmax,ymax]=raw,scale=Math.min((box.width-48)/Math.max(1,xmax-xmin),(box.height-48)/Math.max(1,ymax-ymin));const tx=(box.width-(xmax-xmin)*scale)/2,ty=(box.height-(ymax-ymin)*scale)/2;view={project:(x,y)=>[tx+(x-xmin)*scale,ty+(y-ymin)*scale],scale};const wanted=selected?new Set(selected.selection.cell_ids):null;cellPoints=[];for(const c of units()){const p=view.project(c.x,c.y),active=!wanted||wanted.has(c.cell_id);ctx.globalAlpha=active?1:.12;ctx.fillStyle=colors[c.lineage];ctx.beginPath();if(c.polygon&&el('polygons').checked){c.polygon.forEach((v,i)=>{const q=view.project(v[0],v[1]);if(i)ctx.lineTo(...q);else ctx.moveTo(...q)});ctx.closePath();ctx.fill();ctx.strokeStyle='#ffffff';ctx.lineWidth=.5;ctx.stroke()}else{ctx.arc(p[0],p[1],active?2.6:1.7,0,Math.PI*2);ctx.fill()}cellPoints.push({p,c});}ctx.globalAlpha=1;candidatePoints=[];candidates().forEach((c,i)=>{const p=view.project(...c.center_um);candidatePoints.push({p,c});ctx.beginPath();ctx.arc(p[0],p[1],selected&&c.candidate_id===selected.candidate_id?11:8,0,Math.PI*2);ctx.fillStyle='#fff';ctx.fill();ctx.strokeStyle=selected&&c.candidate_id===selected.candidate_id?'#0c765e':'#263f36';ctx.lineWidth=2;ctx.stroke();ctx.fillStyle='#17372d';ctx.font='10px system-ui';ctx.textAlign='center';ctx.textBaseline='middle';ctx.fillText(String(i+1),p[0],p[1]);});if(selected){const b=selected.selection.window.bounds_um,a=view.project(b[0],b[1]),z=view.project(b[2],b[3]);ctx.setLineDash([5,4]);ctx.strokeStyle='#147b69';ctx.lineWidth=1;ctx.strokeRect(a[0],a[1],z[0]-a[0],z[1]-a[1]);ctx.setLineDash([])}}
canvas.onclick=e=>{const b=canvas.getBoundingClientRect(),x=e.clientX-b.left,y=e.clientY-b.top;const near=candidatePoints.map(p=>({...p,d:Math.hypot(p.p[0]-x,p.p[1]-y)})).sort((a,b)=>a.d-b.d)[0];if(near&&near.d<=16)choose(near.c)};
canvas.onmousemove=e=>{const b=canvas.getBoundingClientRect(),x=e.clientX-b.left,y=e.clientY-b.top;let near=null;for(const point of cellPoints){const d=Math.hypot(point.p[0]-x,point.p[1]-y);if(!near||d<near.d)near={...point,d}};el('pointer').textContent=near&&near.d<10?`${near.c.cell_id} · ${near.c.cell_type} / ${near.c.lineage} · ${near.c.transcript_total} transcripts · QC ${near.c.qc} · ${near.c.geometry}`:''};
download.onclick=()=>{if(!selected)return;const blob=new Blob([JSON.stringify(selected.selection,null,2)+'\n'],{type:'application/json'}),url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download=selected.candidate_id+'.selection.json';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000)};
sample.onchange=list;el('clear').onclick=()=>choose(null);el('polygons').onchange=draw;el('zoom').onchange=draw;window.addEventListener('resize',draw);list();
</script></html>'''
