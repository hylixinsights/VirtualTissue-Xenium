/* The original VirtualTissue dark player, adapted to saved v2/v3 states.
 * Presentation never changes simulation records. Lazy trajectories retain two chunks.
 */
'use strict';
const $=id=>document.getElementById(id);
const data=JSON.parse($('run-data').textContent),frames=data.frames,report=data.report||{};
function inlineMetadata(){
 const lineages=new Set(),order=new Set();let xmin=Infinity,xmax=-Infinity,ymin=Infinity,ymax=-Infinity,polygons=false;
 for(const frame of frames)for(const [id,c] of Object.entries(frame.cells)){
  lineages.add(c.lineage);order.add(id);
  const shape=c.t0_geometry?.status==='valid'?c.t0_geometry:null;
  polygons=polygons||!!shape;
  const contour=shape?shape.polygon_um.map(p=>[p[0]+c.x-shape.anchor_um[0],p[1]+c.y-shape.anchor_um[1]]):[];
  for(const p of [[c.x,c.y],...contour])if(p.every(Number.isFinite)){xmin=Math.min(xmin,p[0]);xmax=Math.max(xmax,p[0]);ymin=Math.min(ymin,p[1]);ymax=Math.max(ymax,p[1])}
 }
 return {lineages:[...lineages].sort(),cell_order:[...order],bounds:Number.isFinite(xmin)?[xmin,xmax,ymin,ymax]:[0,1,0,1],has_polygons:polygons};
}
const meta=data.replay||inlineMetadata(),lazy=!!data.replay;
const slider=$('slider'),canvas=$('canvas'),ctx=canvas.getContext('2d'),play=$('play'),tracks=$('tracks');
const spatial=report.preflight?.interaction_mode!=='well_mixed',lattice=report.kernel==='lattice';
tracks.disabled=!spatial||lattice;
document.title=`VirtualTissue ${report.schema_version?.startsWith('2')?'v2':'v3'} · research replay`;
$('subtitle').textContent=`${report.schema_version?.startsWith('2')?'v2':'v3'} · ${report.mode||'saved'} trajectory · ${lattice?'fixed measured bins':spatial?'spatial tissue':'well-mixed schematic'}`;
$('badge').textContent=`${Number(frames.at(-1)?.time_hours||0).toFixed(2)} h simulated · seed ${report.seed??'—'} · ${report.completed===false?'incomplete run':'recorded run'}`;
slider.max=Math.max(0,frames.length-1);
const esc=s=>String(s??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const human=s=>String(s??'').replaceAll('_',' ');
const PALETTE={malignant:'#ff8a5c',myeloid:'#f2c14e',T_reg:'#cfa2f5',NK:'#38d6cf',stromal:'#b9cf8c',epithelial:'#d8573f',Tumor:'#d8573f',Tumor_proliferating:'#ff8a5c',T_CD8:'#38d6cf',T_CD8_cytotoxic:'#38d6cf',T_CD4:'#6fa8ff',B:'#8f7cf6',B_cell:'#8f7cf6',plasma:'#cfa2f5',Plasma_cell:'#cfa2f5',macrophage:'#f2c14e',Macrophage:'#f2c14e',Macrophage_SPP1:'#f2c14e',fibroblast:'#b9cf8c',Fibroblast:'#b9cf8c',CAF:'#86b35a',endothelial:'#ec6fae',Endothelial:'#ec6fae',unknown:'#5d6684'};
const fallback=['#cfa2f5','#86b35a','#ec6fae','#6fa8ff','#b9cf8c','#f2c14e','#8f7cf6'];
const colors=Object.fromEntries(meta.lineages.map((t,i)=>[t,PALETTE[t]||PALETTE[Object.keys(PALETTE).find(k=>k.toLowerCase()===String(t).toLowerCase())]||fallback[i%fallback.length]]));
const EVENT_STYLE={kill:['CD8 kill','#ff4d5e'],apoptosis:['Apoptosis','#c878ff'],division:['Cell division','#ffab7e'],clear_corpse:['Corpse cleared','#f2c14e'],corpse_degraded:['Corpse removed','#9a8fa6'],deactivate:['Activation ended','#8d97b3'],move:['Movement','#21e6a8'],activate:['Activation','#38d6cf'],exhaust:['Exhaustion','#f2c14e'],recover:['Recovery','#6fa8ff'],population_cap_blocked:['Division blocked by population cap','#ffd36b'],conflict_blocked:['Conflicting action blocked','#8d97b3']};
const warning=lattice?'Fixed measured bins. Occupancy is an abstract amount per bin, not a cell count. Empty bins are faint outlines.':spatial?(meta.has_polygons?'Provided t0 segmentation is shown when valid. Movement rigidly translates outlines; daughter outlines are inherited abstractions. Invalid outlines use centroid markers.':'Centroid markers show recorded positions. Marker size is illustrative. Daughter offsets are synthetic; mechanics and physical contact are not modelled.'):'Schematic layout for well-mixed interactions. Positions are display coordinates and have no anatomical meaning.';
$('warning').textContent=warning;
$('summary').textContent=JSON.stringify({mode:report.mode,completed:report.completed,steps:report.completed_steps,seed:report.seed,events:report.events,usage:report.usage},null,2);
if(lattice){$('tumorLabel').textContent='occupied bins';$('corpseLabel').textContent='empty bins';$('livingLabel').textContent='total bins';$('chartTitle').textContent='Occupied bins over time'}
const index=Object.fromEntries(meta.cell_order.map((id,i)=>[id,i])),nCells=meta.cell_order.length,nColumns=Math.max(1,Math.ceil(Math.sqrt(nCells)));
const [xmin,xmax,ymin,ymax]=spatial?meta.bounds:[0,Math.max(0,nColumns-1),0,Math.max(0,Math.ceil(nCells/nColumns)-1)];
const cache=new Map(),pending=new Map(),focused=new Set();let loadChain=Promise.resolve(),requestVersion=0,currentFrame=null,points=[],selected=null,playing=false,timer=null,drag=null,suppressClick=false;
let view={x:(xmin+xmax)/2,y:(ymin+ymax)/2,scale:1,initialScale:1},viewFitted=false;
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
 document.head.append(script);return promise;
}
function frameAt(position,version){loadChain=loadChain.catch(()=>{}).then(async()=>{if(version!==requestVersion)return null;if(!lazy)return frames[position];const stub=frames[position];const chunk=await loadChunk(stub.chunk);return version===requestVersion?chunk[stub.offset]:null});return loadChain}
function xy(c){const i=index[c.cell_id]??0;return spatial?[c.x,c.y]:[i%nColumns,Math.floor(i/nColumns)]}
function toScreen(p,box){return [(p[0]-view.x)*view.scale+box.width/2,(p[1]-view.y)*view.scale+box.height/2]}
function fitView(){const box=canvas.getBoundingClientRect();view.x=(xmin+xmax)/2;view.y=(ymin+ymax)/2;view.scale=Math.min(Math.max(1,box.width-70)/Math.max(1,xmax-xmin),Math.max(1,box.height-110)/Math.max(1,ymax-ymin));view.initialScale=view.scale;viewFitted=true}
function drawTracks(f,box){
 if(!tracks.checked||!spatial||lattice)return;
 ctx.globalAlpha=.3;ctx.lineWidth=1;
 const chunks=lazy?[...cache.values()]:[frames];
 for(const chunk of chunks)for(const previous of chunk){
  if(previous.step>f.step||previous.step<f.step-100)continue;
  for(const event of previous.events||[]){
   if(event.kind!=='move'||!event.from_xy||!event.to_xy)continue;
   const a=toScreen(event.from_xy,box),b=toScreen(event.to_xy,box),cell=previous.cells[event.cell_id];
   ctx.beginPath();ctx.moveTo(...a);ctx.lineTo(...b);ctx.strokeStyle=colors[cell?.lineage]||'#8d97b3';ctx.stroke();
  }
 }
 ctx.globalAlpha=1;
}
function pointInside(x,y,ring){let inside=false;for(let i=0,j=ring.length-1;i<ring.length;j=i++){const a=ring[i],b=ring[j];if((a[1]>y)!==(b[1]>y)&&x<(b[0]-a[0])*(y-a[1])/(b[1]-a[1])+a[0])inside=!inside}return inside}
function stateColor(c){return c.status==='dead'?'#4a4458':c.exhausted?'#f2c14e':c.activated?'#38d6cf':'#6f7eab'}
function cellColor(c,events){const mode=$('colorMode').value;if(mode==='state')return stateColor(c);if(mode==='event'){const event=events.findLast(e=>e.cell_id===c.cell_id||e.target===c.cell_id);return event?(EVENT_STYLE[event.kind]?.[1]||'#bba6ff'):'#39425f'}if(mode==='signal'){const s=Math.max(0,Number(c.signal)||0);return `hsl(${250-200*s/(1+s)},65%,${35+20*s/(1+s)}%)`}return colors[c.lineage]||'#5d6684'}
function pathCell(polygon,x,y,r){ctx.beginPath();if(polygon){ctx.moveTo(...polygon[0]);for(const vertex of polygon.slice(1))ctx.lineTo(...vertex);ctx.closePath()}else ctx.arc(x,y,r,0,Math.PI*2)}
function drawInteractions(f,box){
 if(!$('showFx').checked)return;
 const byId=new Map(points.map(p=>[p.c.cell_id,p]));
 for(const event of f.events||[]){
  if(!['kill','division','apoptosis','clear_corpse'].includes(event.kind))continue;
  const a=byId.get(event.cell_id),b=byId.get(event.target);
  const color=EVENT_STYLE[event.kind]?.[1]||'#bba6ff';
  if(a&&b&&['kill','clear_corpse'].includes(event.kind)){ctx.beginPath();ctx.moveTo(a.x,a.y);ctx.lineTo(b.x,b.y);ctx.strokeStyle=color;ctx.lineWidth=2;ctx.globalAlpha=.85;ctx.stroke()}
  const p=event.kind==='kill'?b:a;
  if(p&&['kill','division','apoptosis','clear_corpse'].includes(event.kind)){ctx.beginPath();ctx.arc(p.x,p.y,p.r+5,0,Math.PI*2);ctx.strokeStyle=color;ctx.lineWidth=1.5;ctx.globalAlpha=.9;ctx.stroke()}
 }
 ctx.globalAlpha=1;
}
function drawHighlights(f,box){
 if(!$('showHighlights').checked||box.width<300)return;
 const events=(f.events||[]).filter(e=>['kill','division','apoptosis','clear_corpse'].includes(e.kind)).slice(-3),byId=new Map(points.map(p=>[p.c.cell_id,p]));
 events.forEach((event,i)=>{const p=byId.get(event.kind==='kill'?event.target:event.cell_id);if(!p)return;const [title,color]=EVENT_STYLE[event.kind],width=Math.min(175,box.width-30),x=Math.max(12,Math.min(box.width-width-12,p.x+30)),y=100+i*42;if(y>box.height-75)return;ctx.globalAlpha=.95;ctx.strokeStyle=color;ctx.lineWidth=1;ctx.beginPath();ctx.moveTo(x,y+14);ctx.lineTo(p.x,p.y);ctx.stroke();ctx.fillStyle='#0a1225';ctx.fillRect(x,y,width,28);ctx.fillStyle=color;ctx.font='600 11px Inter,system-ui,sans-serif';ctx.fillText(title,x+9,y+18,width-18)});
 ctx.globalAlpha=1;
}
function drawHud(f,box){
 ctx.fillStyle='#8d97b3';ctx.font='11px Inter,system-ui,sans-serif';
 ctx.fillText(spatial?(meta.has_polygons?'Saved geometry · rigid translations':'Recorded centroids · illustrative markers'):'Well-mixed · schematic coordinates',18,box.height-20);
 if(spatial){let um=10;for(const value of [1,5,10,20,50,100,200])if(value*view.scale<=140)um=value;const px=um*view.scale;if(px>=4&&px<box.width-30){ctx.fillStyle='#e7ebf5';ctx.fillRect(18,box.height-45,px,2);ctx.fillText(`${um} µm`,18,box.height-51)}}
}
function metricsAt(position){return frames[position]?.metrics||data.metrics?.[position]||null}
function currentMetrics(f){const m=metricsAt(Number(slider.value));if(m)return m;const cells=Object.values(f.cells),lineages={};for(const c of cells)if(c.status==='alive')lineages[c.lineage]=(lineages[c.lineage]||0)+1;return {living:cells.filter(c=>c.status==='alive').length,dead:cells.filter(c=>c.status==='dead').length,removed:cells.filter(c=>c.status==='removed').length,malignant:lineages.malignant||0,occupied:cells.filter(c=>c.occupancy>0).length,empty:cells.filter(c=>c.occupancy===0).length,event_count:(f.events||[]).length,lineages,events:null}}
function drawChart(){
 const cv=$('chart'),box=cv.getBoundingClientRect(),dpr=window.devicePixelRatio||1,c=cv.getContext('2d');cv.width=box.width*dpr;cv.height=box.height*dpr;c.setTransform(dpr,0,0,dpr,0,0);c.clearRect(0,0,box.width,box.height);
 const key=lattice?'occupied':'malignant',series=frames.map((f,i)=>[f.time_hours,metricsAt(i)?.[key]]),valid=series.filter(s=>Number.isFinite(s[1]));
 if(!valid.length){c.fillStyle='#8d97b3';c.font='11px system-ui';c.fillText('No saved population summary.',8,28);return}
 let max=1;for(const s of valid)max=Math.max(max,s[1]);max*=1.1;
 const end=Math.max(.01,frames.at(-1).time_hours),X=t=>30+(box.width-40)*t/end,Y=n=>box.height-20-(box.height-30)*n/max;
 c.strokeStyle='#232c47';c.lineWidth=1;c.beginPath();c.moveTo(30,8);c.lineTo(30,box.height-20);c.lineTo(box.width-10,box.height-20);c.stroke();
 c.fillStyle='#8d97b3';c.font='10px Inter,system-ui';c.fillText(String(Math.ceil(max)),2,15);c.fillText('0 h',24,box.height-3);c.fillText(`${end.toFixed(1)} h`,Math.max(50,box.width-45),box.height-3);
 c.beginPath();let started=false;for(const s of series){if(!Number.isFinite(s[1])){started=false;continue}if(started)c.lineTo(X(s[0]),Y(s[1]));else{c.moveTo(X(s[0]),Y(s[1]));started=true}}c.strokeStyle=lattice?'#38d6cf':'#ff8a5c';c.lineWidth=2;c.stroke();
 c.beginPath();c.moveTo(X(currentFrame?.time_hours||0),8);c.lineTo(X(currentFrame?.time_hours||0),box.height-20);c.strokeStyle='#e7ebf5';c.lineWidth=1;c.stroke();
}
function drawStrip(){const cv=$('strip'),box=cv.getBoundingClientRect(),dpr=window.devicePixelRatio||1,c=cv.getContext('2d');cv.width=box.width*dpr;cv.height=12*dpr;c.setTransform(dpr,0,0,dpr,0,0);c.clearRect(0,0,box.width,12);const end=Math.max(.01,frames.at(-1).time_hours);let previous={};for(let i=0;i<frames.length;i++){const totals=metricsAt(i)?.events||{};for(const [kind,color,y] of [['division','#ff8a5c',0],['apoptosis','#c878ff',4],['kill','#ff4d5e',4],['clear_corpse','#f2c14e',8]])if((totals[kind]||0)>(previous[kind]||0)){c.fillStyle=color;c.fillRect(8+(box.width-16)*frames[i].time_hours/end,y,1.5,3.5)}previous=totals}}
function bindEventButtons(box,items){for(const button of box.querySelectorAll('button'))button.onclick=()=>{const event=items[Number(button.dataset.index)],id=event.kind==='kill'?event.target:event.cell_id,c=currentFrame.cells[id];if(c){selected=id;const p=xy(c);view.x=p[0];view.y=p[1];draw(currentFrame)}}}
function updateEvents(f){
 const events=f.events||[],items=events.filter(e=>e.kind!=='move').slice(-3).reverse(),highlights=$('highlights');
 highlights.innerHTML=items.length?items.map((e,i)=>`<button class="highlight" data-index="${i}" style="border-left:3px solid ${EVENT_STYLE[e.kind]?.[1]||'#8d97b3'}">${esc(EVENT_STYLE[e.kind]?.[0]||human(e.kind))}<small>cell ${esc(e.cell_id)}${e.target?` → ${esc(e.target)}`:''}${e.daughter?` → daughter ${esc(e.daughter)}`:''} · step ${f.step}</small></button>`).join(''):'<div class="hint">No recorded events at this step. Move along the timeline to explore the trajectory.</div>';
 bindEventButtons(highlights,items);
 const history=(lazy?[...cache.values()].flat():frames).filter(frame=>frame.step<=f.step&&frame.step>=f.step-20).sort((a,b)=>a.step-b.step),recent=[];
 for(const previous of history)for(const event of previous.events||[])if(event.kind!=='move')recent.push({event,time:event.time_hours??previous.time_hours});
 $('feed').innerHTML=recent.length?recent.slice(-12).reverse().map(({event:e,time})=>`<div class="ev"><time>${Number(time).toFixed(2)} h</time><i style="background:${EVENT_STYLE[e.kind]?.[1]||'#8d97b3'}"></i><span>${esc(EVENT_STYLE[e.kind]?.[0]||human(e.kind))} · ${esc(e.cell_id)}${e.target?` → ${esc(e.target)}`:''}${e.daughter?` → ${esc(e.daughter)}`:''}</span></div>`).join(''):'<div class="hint">No events in recently loaded steps.</div>';
}
function updateInspector(f){
 if(selected===null)return;
 const c=f.cells[selected];$('detail').textContent=JSON.stringify(c||{status:'not born'},null,2);
 if(!c){$('inspector').innerHTML='<div class="hint">This cell is absent from the selected saved frame.</div>';return}
 const rows=[['Status',c.status],['Cell type',c.cell_type],['Lineage',human(c.lineage)],['Activated',c.activated],['Exhausted',c.exhausted],['Generation',c.generation],['Parent',c.parent_id??c.shape_inherited_from],['Occupancy',c.occupancy],['Geometry',human(c.geometry_model)],['Activated until',c.activated_until===undefined?undefined:`${Number(c.activated_until).toFixed(2)} h`],['Moving until',c.moving_until===undefined?undefined:`${Number(c.moving_until).toFixed(2)} h`]].filter(([,v])=>v!==undefined&&v!==null);
 const signals=[['Inflammatory signal',c.signal],['Suppressive signal',c.suppressive_signal]].filter(([,v])=>Number.isFinite(v));
 $('inspector').innerHTML=`<div class="title"><i style="background:${colors[c.lineage]}"></i>${esc(c.cell_type||human(c.lineage))}</div><div class="meta">${esc(c.cell_id)} · saved at ${Number(f.time_hours).toFixed(2)} h</div><dl class="kv">${rows.map(([k,v])=>`<dt>${esc(k)}</dt><dd>${esc(String(v))}</dd>`).join('')}${spatial?`<dt>Centroid</dt><dd>${Number(c.x).toFixed(1)}, ${Number(c.y).toFixed(1)} µm</dd>`:''}</dl>${signals.length?'<div class="hint">Saved abstract signal channels (arbitrary units)</div><div class="bars">'+signals.map(([k,v])=>`<div class="bar"><span class="name">${esc(k)}</span><span class="track"><span class="fill" style="display:block;width:${Math.max(0,v)/(1+Math.max(0,v))*100}%;background:#3fd0c9"></span></span><span class="num">${v.toFixed(2)}</span></div>`).join('')+'</div>':''}`;
}
function updatePanel(f){
 const m=currentMetrics(f),totals=m.events||{},value=v=>Number.isFinite(v)?Number(v).toLocaleString('en-US'):'—';
 $('tTumor').textContent=value(lattice?m.occupied:m.malignant);$('tDiv').textContent=value(m.events?totals.division||0:null);$('tDeath').textContent=value(m.events?(totals.kill||0)+(totals.apoptosis||0):null);$('tEat').textContent=value(m.events?totals.clear_corpse||0:null);$('tCorpse').textContent=value(lattice?m.empty:m.dead);$('tLiving').textContent=value(lattice?m.occupied+m.empty:m.living);
 $('deathBreakdown').textContent=m.events?`${totals.kill||0} CD8 kills · ${totals.apoptosis||0} apoptoses · cumulative recorded events`:'';
 $('stats').textContent=lattice?`${m.occupied} occupied bins · ${m.empty} empty bins · ${(f.events||[]).length} events this step`:`${m.living} living · ${m.dead} dead · ${(f.events||[]).length} events this step`;
 const mode=$('colorMode').value;$('legendTitle').textContent=mode==='type'?'Cell types · click to focus':mode==='state'?'Saved cell state':mode==='event'?'Recorded events this step':'Abstract inflammatory signal';
 let entries;if(mode==='state')entries=[['Alive / resting','#6f7eab'],['Activated','#38d6cf'],['Exhausted','#f2c14e'],['Dead','#4a4458']];else if(mode==='event')entries=Object.values(EVENT_STYLE).filter((_,i)=>i<8);else if(mode==='signal')entries=[['Low signal','#5b1f94'],['Higher signal','#57bfe0']];else entries=meta.lineages.map(t=>[human(t),colors[t],t]);
 $('legend').innerHTML=entries.map(([label,color,lineage])=>`<div class="leg ${lineage&&focused.size&&!focused.has(lineage)?'off':''}"${lineage?` data-lineage="${esc(lineage)}"`:''}><i style="background:${color}"></i>${esc(label)}${lineage?`<span class="n">${value(m.lineages?.[lineage]||0)}</span>`:''}</div>`).join('');
 for(const el of $('legend').querySelectorAll('[data-lineage]'))el.onclick=()=>{const lineage=el.dataset.lineage;if(focused.has(lineage))focused.delete(lineage);else focused.add(lineage);draw(f)};
 $('clock').textContent=`${Number(f.time_hours).toFixed(2)} h`;$('time').textContent=`Step ${f.step} · ${Number(f.time_hours).toFixed(2)} simulated hours`;
 $('phase').textContent=`${Number(slider.value)+1} / ${frames.length} saved frames`;
 updateEvents(f);updateInspector(f);drawChart();
}
function draw(f){
 currentFrame=f;if(!viewFitted)fitView();const box=canvas.getBoundingClientRect(),dpr=window.devicePixelRatio||1;canvas.width=box.width*dpr;canvas.height=box.height*dpr;ctx.setTransform(dpr,0,0,dpr,0,0);ctx.clearRect(0,0,box.width,box.height);points=[];
 const gradient=ctx.createRadialGradient(box.width/2,box.height/2,Math.min(box.width,box.height)*.1,box.width/2,box.height/2,Math.max(box.width,box.height)*.75);gradient.addColorStop(0,'#111830');gradient.addColorStop(1,'#070a14');ctx.fillStyle=gradient;ctx.fillRect(0,0,box.width,box.height);drawTracks(f,box);
 const events=f.events||[];
 for(const c of Object.values(f.cells)){
  if(c.status==='removed'&&!lattice)continue;
  const p=xy(c),[x,y]=toScreen(p,box),r=Math.max(1,Math.min(25,(nCells>400?3:5)*view.scale/view.initialScale));
  const shape=spatial&&c.t0_geometry?.status==='valid'?c.t0_geometry:null;
  const polygon=shape?shape.polygon_um.map(p=>toScreen([p[0]+c.x-shape.anchor_um[0],p[1]+c.y-shape.anchor_um[1]],box)):null;
  pathCell(polygon,x,y,r);ctx.fillStyle=c.status==='dead'?'#4a4458':cellColor(c,events);ctx.globalAlpha=(focused.size&&!focused.has(c.lineage)?.18:1)*(lattice&&c.occupancy===0?.15:c.status==='dead'?.45:.9);
  if(lattice&&c.occupancy===0){ctx.strokeStyle=colors[c.lineage];ctx.lineWidth=.8;ctx.stroke()}else ctx.fill();
  ctx.globalAlpha=.55;ctx.strokeStyle=c.status==='dead'?'#9a8fa6':'rgba(255,255,255,.35)';ctx.lineWidth=.65;if(c.status==='dead'||c.geometry_model==='synthetic_inherited_rigid_translation')ctx.setLineDash([3,2]);ctx.stroke();ctx.setLineDash([]);ctx.globalAlpha=1;
  if(c.cell_id===selected){ctx.strokeStyle='#ffffff';ctx.lineWidth=2.4;ctx.stroke()}
  if(!lattice&&c.status==='alive'&&r>2){ctx.beginPath();ctx.arc(x,y,Math.max(1,r*.28),0,Math.PI*2);ctx.fillStyle='#0a0e1a';ctx.globalAlpha=.5;ctx.fill();ctx.globalAlpha=1}
  if(c.status==='dead'){ctx.beginPath();ctx.moveTo(x-r*.6,y-r*.6);ctx.lineTo(x+r*.6,y+r*.6);ctx.strokeStyle='#b8a9c3';ctx.lineWidth=1;ctx.stroke()}
  points.push({x,y,r,c,polygon});
 }
 drawInteractions(f,box);drawHighlights(f,box);drawHud(f,box);updatePanel(f);canvas.dataset.frameStep=String(f.step);document.body.dataset.replayStatus='ready';
}
async function render(){const version=++requestVersion,position=Number(slider.value);document.body.dataset.replayStatus='loading';try{const frame=await frameAt(position,version);if(frame&&version===requestVersion){$('warning').textContent=warning;draw(frame)}}catch(error){if(version===requestVersion){document.body.dataset.replayStatus='error';$('warning').textContent=`Cannot load replay: ${error.message}`;playing=false;clearTimeout(timer);play.textContent='Play'}}}
slider.oninput=render;tracks.onchange=()=>{if(currentFrame)draw(currentFrame)};
function hit(x,y){const containing=points.filter(p=>p.polygon&&pointInside(x,y,p.polygon)),candidates=containing.length?containing:points;let near=null,dist=Infinity;for(const p of candidates){const d=(p.x-x)**2+(p.y-y)**2;if(d<dist){near=p;dist=d}}return near&&(containing.length||dist<Math.max(225,(near.r+5)**2))?near:null}
canvas.onclick=e=>{if(suppressClick){suppressClick=false;return}const b=canvas.getBoundingClientRect(),near=hit(e.clientX-b.left,e.clientY-b.top);if(near){selected=near.c.cell_id;if(currentFrame)draw(currentFrame)}};
async function tick(){if(!playing)return;if(Number(slider.value)>=frames.length-1){playing=false;play.textContent='Play';stopRecording();return}slider.value=Number(slider.value)+1;await render();if(playing)timer=setTimeout(tick,450/Math.max(.1,Number($('speed').value)||1))}
play.onclick=()=>{playing=!playing;play.textContent=playing?'Pause':'Play';if(playing){if(Number(slider.value)>=frames.length-1)slider.value=0;timer=setTimeout(tick,0)}else clearTimeout(timer)};
$('restart').onclick=()=>{playing=false;clearTimeout(timer);play.textContent='Play';slider.value=0;render()};
$('resetView').onclick=()=>{fitView();if(currentFrame)draw(currentFrame)};
let recorder=null;
function stopRecording(){if(recorder&&recorder.state!=='inactive')recorder.stop()}
async function startRecording(){
 if(!window.MediaRecorder||!canvas.captureStream){$('warning').textContent='Video recording is unavailable in this browser. Open this replay in a browser with MediaRecorder support.';return}
 const stream=canvas.captureStream(30),mime=MediaRecorder.isTypeSupported('video/webm;codecs=vp9')?'video/webm;codecs=vp9':'video/webm',chunks=[];
 try{
  recorder=new MediaRecorder(stream,{mimeType:mime,videoBitsPerSecond:8000000});
  recorder.ondataavailable=e=>{if(e.data.size)chunks.push(e.data)};
  recorder.onstop=()=>{const url=URL.createObjectURL(new Blob(chunks,{type:'video/webm'})),link=document.createElement('a');link.href=url;link.download=`virtualtissue_${report.mode||'replay'}_seed_${report.seed??0}.webm`;link.click();setTimeout(()=>URL.revokeObjectURL(url),1000);stream.getTracks().forEach(track=>track.stop());recorder=null;$('record').textContent='● Record video'};
  recorder.start();$('record').textContent='■ Stop recording';playing=false;clearTimeout(timer);slider.value=0;await render();playing=true;play.textContent='Pause';timer=setTimeout(tick,450);
 }catch(error){stream.getTracks().forEach(track=>track.stop());recorder=null;$('warning').textContent=`Cannot record video: ${error.message}`}
}
$('record').onclick=()=>recorder?stopRecording():startRecording();
for(const id of ['colorMode','showFx','showHighlights'])$(id).onchange=()=>{if(currentFrame)draw(currentFrame)};
canvas.addEventListener('wheel',e=>{e.preventDefault();const box=canvas.getBoundingClientRect(),x=e.clientX-box.left,y=e.clientY-box.top,wx=(x-box.width/2)/view.scale+view.x,wy=(y-box.height/2)/view.scale+view.y,factor=Math.exp(-e.deltaY*.0015);view.scale=Math.max(view.initialScale*.2,Math.min(view.initialScale*30,view.scale*factor));view.x=wx-(x-box.width/2)/view.scale;view.y=wy-(y-box.height/2)/view.scale;if(currentFrame)draw(currentFrame)},{passive:false});
canvas.addEventListener('mousedown',e=>{drag={x:e.clientX,y:e.clientY,vx:view.x,vy:view.y,moved:false};canvas.classList.add('dragging')});
window.addEventListener('mouseup',()=>{if(drag?.moved)suppressClick=true;drag=null;canvas.classList.remove('dragging')});
canvas.addEventListener('mousemove',e=>{const box=canvas.getBoundingClientRect(),x=e.clientX-box.left,y=e.clientY-box.top;if(drag){const dx=e.clientX-drag.x,dy=e.clientY-drag.y;if(Math.hypot(dx,dy)>3)drag.moved=true;view.x=drag.vx-dx/view.scale;view.y=drag.vy-dy/view.scale;if(currentFrame)draw(currentFrame);return}const p=hit(x,y),tip=$('tooltip');if(!p){tip.style.display='none';return}tip.innerHTML=`<b>${esc(p.c.cell_type||human(p.c.lineage))}</b><br>${esc(p.c.cell_id)} · ${esc(p.c.status)}${p.c.activated?' · activated':''}${p.c.exhausted?' · exhausted':''}`;tip.style.display='block';tip.style.left=Math.max(8,Math.min(x+14,box.width-220))+'px';tip.style.top=Math.max(8,Math.min(y+14,box.height-70))+'px'});
canvas.addEventListener('mouseleave',()=>{$('tooltip').style.display='none'});
canvas.addEventListener('dblclick',()=>{fitView();if(currentFrame)draw(currentFrame)});
window.addEventListener('resize',()=>{drawStrip();if(currentFrame)draw(currentFrame)});
window.addEventListener('keydown',e=>{if(e.code==='Space'&&!['INPUT','SELECT','BUTTON','SUMMARY'].includes(e.target.tagName)){e.preventDefault();play.onclick()}});
drawStrip();cacheDiagnostics();render();
