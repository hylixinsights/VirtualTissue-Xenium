"""Dependency-free offline replay. Well-mixed layouts are explicitly schematic."""
from pathlib import Path
import json


def write_replay(path, frames, report):
    data = json.dumps(dict(frames=frames, report=report), allow_nan=False).replace("<", "\\u003c")
    template = r'''<!doctype html><html lang="en"><meta charset="utf-8">
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
slider.max=frames.length-1;const types=[...new Set(frames.flatMap(f=>Object.values(f.cells).map(c=>c.lineage)))].sort();
const colors=Object.fromEntries(types.map((t,i)=>[t,`hsl(${(i*137.5)%360},48%,44%)`]));
for(const t of types){const label=document.createElement('span'),dot=document.createElement('span');dot.className='dot';dot.style.background=colors[t];label.append(dot,document.createTextNode(t));document.getElementById('legend').append(label)}
const spatial=report.preflight.interaction_mode==='spatial';
document.getElementById('warning').textContent=spatial?'Centroid graph. Daughter offsets are synthetic; mechanics and physical contact are not modelled.':'Schematic layout for well-mixed interactions. Positions are display coordinates and have no anatomical meaning.';
document.getElementById('summary').textContent=JSON.stringify({mode:report.mode,completed:report.completed,steps:report.completed_steps,seed:report.seed,events:report.events,usage:report.usage},null,2);
const all=Object.values(frames.at(-1).cells),index=Object.fromEntries(all.map((c,i)=>[c.cell_id,i]));
const positions=all.map(c=>spatial?[c.x,c.y]:[index[c.cell_id]%Math.ceil(Math.sqrt(all.length)),Math.floor(index[c.cell_id]/Math.ceil(Math.sqrt(all.length)))]);
const xs=positions.map(p=>p[0]),ys=positions.map(p=>p[1]),xmin=Math.min(...xs),xmax=Math.max(...xs),ymin=Math.min(...ys),ymax=Math.max(...ys);
let points=[],selected=null,timer=null;
function render(){const f=frames[Number(slider.value)],box=canvas.getBoundingClientRect(),dpr=devicePixelRatio||1;canvas.width=box.width*dpr;canvas.height=box.height*dpr;ctx.setTransform(dpr,0,0,dpr,0,0);ctx.clearRect(0,0,box.width,box.height);points=[];
document.getElementById('time').textContent=`Step ${f.step} · ${f.time_hours.toFixed(2)} assumed hours`;
const cells=Object.values(f.cells);document.getElementById('stats').textContent=`${cells.filter(c=>c.status==='alive').length} living · ${cells.filter(c=>c.status==='dead').length} dead · ${f.events.length} events this step`;
for(const c of cells){if(c.status==='removed')continue;const p=positions[index[c.cell_id]],x=25+(p[0]-xmin)/Math.max(1,xmax-xmin)*(box.width-50),y=25+(p[1]-ymin)/Math.max(1,ymax-ymin)*(box.height-50),r=all.length>400?3:5;ctx.beginPath();ctx.arc(x,y,r,0,Math.PI*2);ctx.fillStyle=colors[c.lineage];ctx.globalAlpha=c.status==='dead'?.35:1;ctx.fill();ctx.globalAlpha=1;if(c.cell_id===selected){ctx.strokeStyle='#111';ctx.lineWidth=2;ctx.stroke()}if(c.status==='dead'){ctx.beginPath();ctx.moveTo(x-r,y-r);ctx.lineTo(x+r,y+r);ctx.strokeStyle='#222';ctx.lineWidth=1;ctx.stroke()}points.push({x,y,c})}
if(selected)document.getElementById('detail').textContent=JSON.stringify(f.cells[selected]||{status:'not born'},null,2)}
slider.oninput=render;canvas.onclick=e=>{const b=canvas.getBoundingClientRect(),x=e.clientX-b.left,y=e.clientY-b.top;const near=points.map(p=>({...p,d:(p.x-x)**2+(p.y-y)**2})).sort((a,b)=>a.d-b.d)[0];if(near&&near.d<225){selected=near.c.cell_id;render()}};
document.getElementById('play').onclick=()=>{if(timer){clearInterval(timer);timer=null;document.getElementById('play').textContent='Play';return}document.getElementById('play').textContent='Pause';timer=setInterval(()=>{slider.value=(Number(slider.value)+1)%frames.length;render()},450)};window.addEventListener('resize',render);render();
</script></html>'''
    Path(path).write_text(template.replace("__DATA__", data), encoding="utf-8")
