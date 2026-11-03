"""Dependency-free offline replay with a streaming writer for v3 trajectories."""
from pathlib import Path
import json
import hashlib


TEMPLATE = r'''<!doctype html><html lang="en"><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>VirtualTissue v2 · research replay</title>
<style>
*{box-sizing:border-box}body{font:15px system-ui;margin:0;background:#f3f6f5;color:#153a38}header{padding:22px 30px;background:#143d3a;color:#fff}
h1{font-size:25px;margin:0 0 8px}header p{margin:0;color:#c9ddda}main{padding:24px;display:grid;grid-template-columns:minmax(400px,1fr) 340px;gap:22px}
.panel{background:#fff;border:1px solid #d3e0dd;border-radius:12px;padding:20px}canvas{width:100%;height:560px;background:#f7faf9;border:1px solid #d3e0dd;border-radius:8px}
input{width:100%;accent-color:#147d77}button{padding:8px 18px;border:0;border-radius:6px;background:#147d77;color:white;cursor:pointer}pre{white-space:pre-wrap;font-size:12px;max-height:480px;overflow:auto}
.muted{color:#657e79;font-size:13px}.legend{display:flex;gap:12px;flex-wrap:wrap;margin:12px 0}.dot{display:inline-block;width:10px;height:10px;border-radius:50%;margin-right:4px}
#stats{font-weight:600;margin:14px 0}#warning{background:#f5edcc;padding:12px;border-radius:6px;margin-bottom:14px;font-size:13px}
@media(max-width:850px){main{display:block;padding:10px}.panel{margin-bottom:12px}canvas{height:420px}}
</style><header><h1>VirtualTissue v2</h1><p>Context-conditioned cell states · inspectable decisions · research preview</p></header>
<main><section class="panel"><div id="warning"></div><button id="play">Play</button> <span id="time"></span>
<input id="slider" type="range" min="0" value="0"><div id="stats"></div><canvas id="canvas"></canvas>
<div class="legend" id="legend"></div><p class="muted">Click a cell to inspect its simulated state. Dead cells are crossed out; cleared cells are hidden.</p></section>
<aside class="panel"><h2>Run evidence</h2><p class="muted">Heuristic preferences are uncalibrated. Expression remains fixed at t0. This replay does not establish treatment efficacy.</p><pre id="summary"></pre><h2>Selected cell</h2><pre id="detail">Click a cell in the replay.</pre></aside></main>
<script id="run-data" type="application/json">__DATA__</script><script>
const data=JSON.parse(document.getElementById('run-data').textContent),frames=data.frames,report=data.report;
const slider=document.getElementById('slider'),canvas=document.getElementById('canvas'),ctx=canvas.getContext('2d');
slider.max=frames.length-1;const typeSet=new Set();for(const f of frames)for(const c of Object.values(f.cells))typeSet.add(c.lineage);const types=[...typeSet].sort();
const colors=Object.fromEntries(types.map((t,i)=>[t,`hsl(${(i*137.5)%360},48%,44%)`]));
for(const t of types){const label=document.createElement('span'),dot=document.createElement('span');dot.className='dot';dot.style.background=colors[t];label.append(dot,document.createTextNode(t));document.getElementById('legend').append(label)}
const spatial=report.preflight.interaction_mode==='spatial',lattice=report.kernel==='lattice';if(report.schema_version==='3.0'){document.title='VirtualTissue v3 · research replay';document.querySelector('h1').textContent='VirtualTissue v3'}
document.getElementById('warning').textContent=lattice?'Fixed measured bins. Occupancy is an abstract amount per bin, not a cell count. Empty bins are faint outlines.':spatial?'Centroid graph. Daughter offsets are synthetic; mechanics and physical contact are not modelled.':'Schematic layout for well-mixed interactions. Positions are display coordinates and have no anatomical meaning.';
document.getElementById('summary').textContent=JSON.stringify({mode:report.mode,completed:report.completed,steps:report.completed_steps,seed:report.seed,events:report.events,usage:report.usage},null,2);
const all=Object.values(frames.at(-1).cells),index=Object.fromEntries(all.map((c,i)=>[c.cell_id,i]));
const positions=all.map(c=>spatial?[c.x,c.y]:[index[c.cell_id]%Math.ceil(Math.sqrt(all.length)),Math.floor(index[c.cell_id]/Math.ceil(Math.sqrt(all.length)))]);
let xmin=Infinity,xmax=-Infinity,ymin=Infinity,ymax=-Infinity;function bound(x,y){xmin=Math.min(xmin,x);xmax=Math.max(xmax,x);ymin=Math.min(ymin,y);ymax=Math.max(ymax,y)}if(spatial){for(const f of frames)for(const c of Object.values(f.cells))bound(c.x,c.y)}else{for(const p of positions)bound(p[0],p[1])}
let points=[],selected=null,timer=null;
function render(){const f=frames[Number(slider.value)],box=canvas.getBoundingClientRect(),dpr=devicePixelRatio||1;canvas.width=box.width*dpr;canvas.height=box.height*dpr;ctx.setTransform(dpr,0,0,dpr,0,0);ctx.clearRect(0,0,box.width,box.height);points=[];
document.getElementById('time').textContent=`Step ${f.step} · ${f.time_hours.toFixed(2)} assumed hours`;
const cells=Object.values(f.cells);document.getElementById('stats').textContent=lattice?`${cells.filter(c=>c.occupancy>0).length} occupied bins · ${cells.filter(c=>c.occupancy===0).length} empty bins · ${f.events.length} events this step`:`${cells.filter(c=>c.status==='alive').length} living · ${cells.filter(c=>c.status==='dead').length} dead · ${f.events.length} events this step`;
for(const c of cells){if(c.status==='removed'&&!lattice)continue;const p=spatial?[c.x,c.y]:positions[index[c.cell_id]],x=25+(p[0]-xmin)/Math.max(1,xmax-xmin)*(box.width-50),y=25+(p[1]-ymin)/Math.max(1,ymax-ymin)*(box.height-50),r=all.length>400?3:5;ctx.beginPath();ctx.arc(x,y,r,0,Math.PI*2);ctx.fillStyle=colors[c.lineage];ctx.globalAlpha=lattice&&c.occupancy===0?.15:c.status==='dead'?.35:1;if(lattice&&c.occupancy===0){ctx.strokeStyle=colors[c.lineage];ctx.lineWidth=1;ctx.stroke()}else{ctx.fill()}ctx.globalAlpha=1;if(c.cell_id===selected){ctx.strokeStyle='#111';ctx.lineWidth=2;ctx.stroke()}if(c.status==='dead'){ctx.beginPath();ctx.moveTo(x-r,y-r);ctx.lineTo(x+r,y+r);ctx.strokeStyle='#222';ctx.lineWidth=1;ctx.stroke()}points.push({x,y,c})}
if(selected)document.getElementById('detail').textContent=JSON.stringify(f.cells[selected]||{status:'not born'},null,2)}
slider.oninput=render;canvas.onclick=e=>{const b=canvas.getBoundingClientRect(),x=e.clientX-b.left,y=e.clientY-b.top;const near=points.map(p=>({...p,d:(p.x-x)**2+(p.y-y)**2})).sort((a,b)=>a.d-b.d)[0];if(near&&near.d<225){selected=near.c.cell_id;render()}};
document.getElementById('play').onclick=()=>{if(timer){clearInterval(timer);timer=null;document.getElementById('play').textContent='Play';return}document.getElementById('play').textContent='Pause';timer=setInterval(()=>{slider.value=(Number(slider.value)+1)%frames.length;render()},450)};window.addEventListener('resize',render);render();
</script></html>'''


def write_replay(path, frames, report):
    data = json.dumps(dict(frames=frames, report=report), allow_nan=False).replace("<", "\\u003c")
    Path(path).write_text(TEMPLATE.replace("__DATA__", data), encoding="utf8")



LAZY_SCRIPT = r'''
const data=JSON.parse(document.getElementById('run-data').textContent),frames=data.frames,report=data.report,meta=data.replay;
const slider=document.getElementById('slider'),canvas=document.getElementById('canvas'),ctx=canvas.getContext('2d'),play=document.getElementById('play');
const spatial=report.preflight.interaction_mode==='spatial',lattice=report.kernel==='lattice',tracks=document.getElementById('tracks');
tracks.disabled=!spatial||lattice;
document.title='VirtualTissue v3 · research replay';document.querySelector('h1').textContent='VirtualTissue v3';
slider.max=frames.length-1;
const types=meta.lineages,colors=Object.fromEntries(types.map((t,i)=>[t,`hsl(${(i*137.5)%360},48%,44%)`]));
for(const t of types){const label=document.createElement('span'),dot=document.createElement('span');dot.className='dot';dot.style.background=colors[t];label.append(dot,document.createTextNode(t));document.getElementById('legend').append(label)}
const warning=lattice?'Fixed measured bins. Occupancy is an abstract amount per bin, not a cell count. Empty bins are faint outlines.':spatial?(meta.has_polygons?'Provided t0 segmentation is shown when valid. Movement rigidly translates outlines; daughter outlines are inherited abstractions. Invalid outlines use centroids.':'Centroid graph. Daughter offsets are synthetic; mechanics and physical contact are not modelled.'):'Schematic layout for well-mixed interactions. Positions are display coordinates and have no anatomical meaning.';
document.getElementById('warning').textContent=warning;
document.getElementById('summary').textContent=JSON.stringify({mode:report.mode,completed:report.completed,steps:report.completed_steps,seed:report.seed,events:report.events,usage:report.usage},null,2);
const index=Object.fromEntries(meta.cell_order.map((id,i)=>[id,i])),nCells=meta.cell_order.length,nColumns=Math.ceil(Math.sqrt(nCells));
const [xmin,xmax,ymin,ymax]=spatial?meta.bounds:[0,Math.max(0,nColumns-1),0,Math.max(0,Math.ceil(nCells/nColumns)-1)];
const cache=new Map(),pending=new Map();let loadChain=Promise.resolve(),requestVersion=0,currentFrame=null,points=[],selected=null,playing=false,timer=null;
function cacheDiagnostics(){document.body.dataset.replayCachedChunks=String(cache.size);document.body.dataset.replayCachedFrames=String([...cache.values()].reduce((n,chunk)=>n+chunk.length,0));document.body.dataset.replayLoadingChunks=String(pending.size)}
function retainChunk(id,chunk){cache.delete(id);cache.set(id,chunk);while(cache.size>2)cache.delete(cache.keys().next().value);cacheDiagnostics()}
window.__vtxReplayChunk=(id,chunk,shapes={})=>{const job=pending.get(id);if(!job)return;if(!Array.isArray(chunk)||chunk.length!==meta.chunks[id].count){job.fail(new Error('Replay chunk length differs from its manifest.'));return}
 for(const frame of chunk)for(const cell of Object.values(frame.cells)){if(cell.t0_geometry_ref){const shape=shapes[cell.t0_geometry_ref];if(!shape){job.fail(new Error('Replay shape reference is missing.'));return}cell.t0_geometry=shape;delete cell.t0_geometry_ref}}
 retainChunk(id,chunk);job.finish(chunk)};
function loadChunk(id){
 if(cache.has(id)){const chunk=cache.get(id);retainChunk(id,chunk);return Promise.resolve(chunk)}
 if(pending.has(id))return pending.get(id).promise;
 const spec=meta.chunks[id];if(!spec)return Promise.reject(new Error('Replay chunk is absent from its manifest.'));
 const script=document.createElement('script');script.src=spec.src;script.async=true;
 let resolve,reject;const promise=new Promise((a,b)=>{resolve=a;reject=b});
 const timeout=setTimeout(()=>job.fail(new Error('Timed out loading the local replay chunk.')),15000);
 const clean=()=>{clearTimeout(timeout);pending.delete(id);script.remove();cacheDiagnostics()};
 const job={promise,finish:chunk=>{clean();resolve(chunk)},fail:error=>{clean();reject(error)}};
 pending.set(id,job);cacheDiagnostics();script.onerror=()=>job.fail(new Error('Keep replay.html together with its replay_frames folder.'));
 script.onload=()=>{if(pending.has(id))job.fail(new Error('The local replay chunk did not supply its data.'))};
 document.head.append(script);return promise
}
function frameAt(position,version){loadChain=loadChain.catch(()=>{}).then(async()=>{if(version!==requestVersion)return null;const stub=frames[position];const chunk=await loadChunk(stub.chunk);return version===requestVersion?chunk[stub.offset]:null});return loadChain}
function drawTracks(f,box){
 if(!tracks.checked||!spatial||lattice)return;
 const xy=p=>[25+(p[0]-xmin)/Math.max(1,xmax-xmin)*(box.width-50),25+(p[1]-ymin)/Math.max(1,ymax-ymin)*(box.height-50)];
 ctx.globalAlpha=.3;ctx.lineWidth=1;
 for(const chunk of cache.values())for(const previous of chunk){
  if(previous.step>f.step||previous.step<f.step-100)continue;
  for(const event of previous.events){
   if(event.kind!=='move'||!event.from_xy||!event.to_xy)continue;
   const a=xy(event.from_xy),b=xy(event.to_xy),cell=previous.cells[event.cell_id];
   ctx.beginPath();ctx.moveTo(...a);ctx.lineTo(...b);ctx.strokeStyle=colors[cell?cell.lineage:'']||'#657e79';ctx.stroke()
  }
 }
 ctx.globalAlpha=1
}
function pointInside(x,y,ring){let inside=false;for(let i=0,j=ring.length-1;i<ring.length;j=i++){const a=ring[i],b=ring[j];if((a[1]>y)!==(b[1]>y)&&x<(b[0]-a[0])*(y-a[1])/(b[1]-a[1])+a[0])inside=!inside}return inside}
function draw(f){
 currentFrame=f;const box=canvas.getBoundingClientRect(),dpr=devicePixelRatio||1;canvas.width=box.width*dpr;canvas.height=box.height*dpr;ctx.setTransform(dpr,0,0,dpr,0,0);ctx.clearRect(0,0,box.width,box.height);points=[];drawTracks(f,box);
 document.getElementById('time').textContent=`Step ${f.step} · ${f.time_hours.toFixed(2)} assumed hours`;
 canvas.dataset.frameStep=String(f.step);document.body.dataset.replayStatus='ready';
 const cells=Object.values(f.cells);document.getElementById('stats').textContent=lattice?`${cells.filter(c=>c.occupancy>0).length} occupied bins · ${cells.filter(c=>c.occupancy===0).length} empty bins · ${f.events.length} events this step`:`${cells.filter(c=>c.status==='alive').length} living · ${cells.filter(c=>c.status==='dead').length} dead · ${f.events.length} events this step`;
 for(const c of cells){
  if(c.status==='removed'&&!lattice)continue;
  const i=index[c.cell_id],p=spatial?[c.x,c.y]:[i%nColumns,Math.floor(i/nColumns)],x=25+(p[0]-xmin)/Math.max(1,xmax-xmin)*(box.width-50),y=25+(p[1]-ymin)/Math.max(1,ymax-ymin)*(box.height-50),r=nCells>400?3:5;
  const shape=spatial&&c.t0_geometry&&c.t0_geometry.status==='valid'?c.t0_geometry:null;
  const polygon=shape?shape.polygon_um.map(p=>[25+(p[0]+c.x-shape.anchor_um[0]-xmin)/Math.max(1,xmax-xmin)*(box.width-50),25+(p[1]+c.y-shape.anchor_um[1]-ymin)/Math.max(1,ymax-ymin)*(box.height-50)]):null;
  ctx.beginPath();if(polygon){ctx.moveTo(...polygon[0]);for(const vertex of polygon.slice(1))ctx.lineTo(...vertex);ctx.closePath()}else ctx.arc(x,y,r,0,Math.PI*2);
  ctx.fillStyle=colors[c.lineage];ctx.globalAlpha=lattice&&c.occupancy===0?.15:c.status==='dead'?.35:polygon?.6:1;
  if(lattice&&c.occupancy===0){ctx.strokeStyle=colors[c.lineage];ctx.lineWidth=1;ctx.stroke()}else ctx.fill();ctx.globalAlpha=1;
  if(polygon){ctx.strokeStyle=colors[c.lineage];ctx.lineWidth=.7;if(c.geometry_model==='synthetic_inherited_rigid_translation')ctx.setLineDash([3,2]);ctx.stroke();ctx.setLineDash([])}
  if(c.cell_id===selected){ctx.strokeStyle='#111';ctx.lineWidth=2;ctx.stroke()}
  if(c.status==='dead'){ctx.beginPath();ctx.moveTo(x-r,y-r);ctx.lineTo(x+r,y+r);ctx.strokeStyle='#222';ctx.lineWidth=1;ctx.stroke()}
  points.push({x,y,c,polygon})
 }
 if(selected)document.getElementById('detail').textContent=JSON.stringify(f.cells[selected]||{status:'not born'},null,2)
}
async function render(){
 const version=++requestVersion,position=Number(slider.value);document.body.dataset.replayStatus='loading';
 try{const frame=await frameAt(position,version);if(frame&&version===requestVersion){document.getElementById('warning').textContent=warning;draw(frame)}}
 catch(error){if(version===requestVersion){document.body.dataset.replayStatus='error';document.getElementById('warning').textContent=`Cannot load replay: ${error.message}`;playing=false;clearTimeout(timer);play.textContent='Play'}}
}
slider.oninput=render;tracks.onchange=()=>{if(currentFrame)draw(currentFrame)};
canvas.onclick=e=>{const b=canvas.getBoundingClientRect(),x=e.clientX-b.left,y=e.clientY-b.top;const containing=points.filter(p=>p.polygon&&pointInside(x,y,p.polygon));const near=(containing.length?containing:points).map(p=>({...p,d:(p.x-x)**2+(p.y-y)**2})).sort((a,b)=>a.d-b.d)[0];if(near&&(containing.length||near.d<225)){selected=near.c.cell_id;if(currentFrame)draw(currentFrame)}};
async function tick(){if(!playing)return;slider.value=(Number(slider.value)+1)%frames.length;await render();if(playing)timer=setTimeout(tick,450)}
play.onclick=()=>{playing=!playing;play.textContent=playing?'Pause':'Play';if(playing)timer=setTimeout(tick,0);else clearTimeout(timer)};
window.addEventListener('resize',()=>{if(currentFrame)draw(currentFrame)});cacheDiagnostics();render();
'''

LAZY_TEMPLATE = (TEMPLATE.split("</script><script>", 1)[0]
    .replace('<span id="time"></span>', '<span id="time"></span> <label class="muted" title="Shows movement recorded in recently viewed steps."><input id="tracks" type="checkbox" style="width:auto;margin-left:12px"> Recent tracks</label>')
    + "</script><script>\n" + LAZY_SCRIPT + "\n</script></html>")


def write_replay_stream(path, frames_jsonl, report, *, chunk_size=50):
    """Write bounded offline chunks; immutable t0 shapes are shared within a chunk.

    The player rehydrates geometry references without copying the shared shape.
    At most two chunks are retained. Full JSONL/JSON snapshots remain the audit
    exports; the sidecars are a lossless presentation representation.
    """
    from ..v3.geometry import projected_polygon
    if type(chunk_size) is not int or not 1 <= chunk_size <= 500:
        raise ValueError("Replay chunk_size must be an integer from 1 to 500.")
    path = Path(path)
    sidecars = path.with_name(path.stem + "_frames")
    sidecars.mkdir(parents=True, exist_ok=True)
    stubs, chunks, lineages, cell_order, seen_ids = [], [], set(), [], set()
    bounds = [None, None, None, None]
    stream, shapes, has_polygons = None, {}, False
    def close_chunk():
        if stream is not None and not stream.closed:
            stream.write("]")
            if shapes:
                stream.write("," + json.dumps(shapes,allow_nan=False).replace("<","\\u003c"))
            stream.write(");\n")
            stream.close()
    try:
        with Path(frames_jsonl).open(encoding="utf8") as source:
            for raw in source:
                if not raw.strip():
                    continue
                frame = json.loads(raw)
                index = len(stubs)
                chunk_id, offset = divmod(index, chunk_size)
                if offset == 0:
                    close_chunk()
                    shapes = {}
                    filename = f"chunk_{chunk_id:05d}.js"
                    stream = (sidecars / filename).open("w", encoding="utf8")
                    stream.write(f"window.__vtxReplayChunk({chunk_id},[\n")
                    chunks.append(dict(src=sidecars.name + "/" + filename, start=index, count=0))
                displayed = dict(frame,cells={})
                for cid, cell in frame["cells"].items():
                    lineages.add(cell["lineage"])
                    if cid not in seen_ids:
                        seen_ids.add(cid)
                        cell_order.append(cid)
                    contour = projected_polygon(cell)
                    has_polygons = has_polygons or contour is not None
                    positions = (contour or []) + ([[cell["x"],cell["y"]]] if "x" in cell and "y" in cell else [])
                    for x,y in positions:
                        bounds[0] = x if bounds[0] is None else min(bounds[0], x)
                        bounds[1] = x if bounds[1] is None else max(bounds[1], x)
                        bounds[2] = y if bounds[2] is None else min(bounds[2], y)
                        bounds[3] = y if bounds[3] is None else max(bounds[3], y)
                    value = dict(cell)
                    if "t0_geometry" in value:
                        geometry = value.pop("t0_geometry")
                        geometry_key = hashlib.sha256(json.dumps(geometry,sort_keys=True,separators=(",",":"),allow_nan=False).encode()).hexdigest()
                        shapes[geometry_key] = geometry
                        value["t0_geometry_ref"] = geometry_key
                    displayed["cells"][cid] = value
                if offset:
                    stream.write(",\n")
                stream.write(json.dumps(displayed, allow_nan=False).replace("<", "\\u003c"))
                chunks[-1]["count"] += 1
                stubs.append(dict(step=frame["step"], time_hours=frame["time_hours"], chunk=chunk_id, offset=offset))
    finally:
        close_chunk()
    if not stubs:
        raise ValueError("Replay requires at least one committed frame.")
    metadata = dict(format="vtx-offline-chunks-2", chunk_size=chunk_size, cache_chunks=2,
                    chunks=chunks, lineages=sorted(lineages), cell_order=cell_order, has_polygons=has_polygons,
                    geometry_encoding="immutable t0 records shared once per chunk; runtime contours are rigid translations",
                    bounds=bounds if bounds[0] is not None else [0, 1, 0, 1])
    value = dict(frames=stubs, report=report, replay=metadata)
    encoded = json.dumps(value, allow_nan=False).replace("<", "\\u003c")
    path.write_text(LAZY_TEMPLATE.replace("__DATA__", encoded), encoding="utf8")
    (sidecars / "manifest.json").write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf8")
