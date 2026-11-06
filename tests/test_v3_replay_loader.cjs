// A DOM test double executes the actual lazy player and its generated sidecars.
// This tests loader/cache/seek behavior, not browser visual rendering.
const assert=require('node:assert/strict'),fs=require('node:fs'),path=require('node:path'),vm=require('node:vm');
const folder=process.argv[2],html=fs.readFileSync(path.join(folder,'replay.html'),'utf8');
const encoded=html.match(/<script id="run-data" type="application\/json">([\s\S]*?)<\/script>/)[1];
const data=JSON.parse(encoded);
const player=html.split('</script><script>')[1].split('</script>')[0];
let points=[],strokes=0,polygonDraws=0,currentPath=[],mapFills=[];const requested=[],chartSeriesLengths=[],polygonVertexCounts={};
function context2d(id){
 let drawing=[],arcs=[];
 return {
  setTransform(){},clearRect(){if(id==='canvas'){points=[];mapFills=[]}},
  beginPath(){drawing=[];arcs=[];if(id==='canvas')currentPath=[]},
  arc(x,y,r){arcs.push({x,y,r});if(id==='canvas')points.push({x,y,r})},
  fill(){if(id==='canvas')mapFills.push({color:this.fillStyle,alpha:this.globalAlpha,vertices:drawing.map(p=>p.slice()),arcs:arcs.map(p=>({...p}))})},fillRect(){},fillText(){},
  stroke(){if(id==='canvas')strokes++;if(id==='chart')chartSeriesLengths.push(drawing.length)},
  moveTo(x,y){drawing=[[x,y]];if(id==='canvas')currentPath=[[x,y]]},
  lineTo(x,y){drawing.push([x,y]);if(id==='canvas')currentPath.push([x,y])},
  closePath(){if(id==='canvas'&&currentPath.length){polygonDraws++;polygonVertexCounts[currentPath.length]=(polygonVertexCounts[currentPath.length]||0)+1;points.push({x:currentPath.reduce((s,p)=>s+p[0],0)/currentPath.length,y:currentPath.reduce((s,p)=>s+p[1],0)/currentPath.length,r:5})}},
  setLineDash(){},createRadialGradient(){return {addColorStop(){}}},createLinearGradient(){return {addColorStop(){}}},
  measureText(t){return {width:String(t).length*6}},roundRect(){},save(){},restore(){},translate(){},scale(){}
 };
}
function element(id=''){const ctx=context2d(id);return {id,dataset:{},style:{},textContent:'',innerHTML:'',value:'0',checked:false,disabled:false,children:[],append(...items){this.children.push(...items)},appendChild(item){this.append(item)},remove(){},setAttribute(k,v){this[k]=v},getAttribute(k){return this[k]},addEventListener(k,fn){this['on'+k]=fn},classList:{add(){},remove(){},toggle(){}},querySelectorAll(){return this.children.filter(x=>x.tagName==='BUTTON')},getBoundingClientRect(){return {width:600,height:id==='chart'?110:560,left:0,top:0,bottom:560}},getContext(){return ctx}}}
const ids={};
for(const tag of html.matchAll(/<[^>]+\bid="([^"]+)"[^>]*>/g)){
 const e=element(tag[1]);e.checked=/\bchecked\b/.test(tag[0]);const value=tag[0].match(/\bvalue="([^"]*)"/);if(value)e.value=value[1];ids[e.id]=e;
}
for(const select of html.matchAll(/<select\b[^>]*\bid="([^"]+)"[^>]*>([\s\S]*?)<\/select>/g)){
 const options=[...select[2].matchAll(/<option\b[^>]*\bvalue="([^"]+)"([^>]*)>/g)];
 const choice=options.find(option=>/\bselected\b/.test(option[2]))||options[0];if(choice)ids[select[1]].value=choice[1];
}
ids['run-data'].textContent=encoded;const body=element('body'),heading=element('h1');
let realm;
const document={
 body,title:'',getElementById:id=>ids[id],querySelector:selector=>heading,
 createElement:tag=>Object.assign(element(tag),{tagName:tag.toUpperCase()}),createTextNode:text=>({textContent:text}),addEventListener(){},
 head:{append(script){
  requested.push(script.src);
  setImmediate(()=>{
   try{vm.runInContext(fs.readFileSync(path.join(folder,script.src),'utf8'),realm);if(script.onload)script.onload()}
   catch(error){if(script.onerror)script.onerror(error)}
  });
 }}
};
realm=vm.createContext({document,devicePixelRatio:1,setTimeout,clearTimeout,console,performance:{now:()=>Date.now()},requestAnimationFrame:fn=>setImmediate(()=>fn(Date.now()))});
realm.window=realm;realm.window.addEventListener=()=>{};
vm.runInContext(player,realm);
async function waitFor(predicate){for(let i=0;i<100;i++){if(predicate())return;await new Promise(resolve=>setTimeout(resolve,5))}throw new Error('Timed out waiting for replay state: '+JSON.stringify(body.dataset))}
async function seek(step){ids.slider.value=String(step);await ids.slider.oninput();assert.equal(ids.canvas.dataset.frameStep,String(step));assert.equal(body.dataset.replayStatus,'ready');assert(Number(body.dataset.replayCachedChunks)<=2);assert(Number(body.dataset.replayCachedFrames)<=100);const metrics=data.frames[step].metrics;if(metrics){assert.equal(Number(ids.tTumor.textContent),metrics.malignant,'Sidebar tumor count follows selected frame');assert.equal(Number(ids.tLiving.textContent),metrics.living,'Sidebar living count follows selected frame');assert.equal(Number(ids.tDiv.textContent),metrics.events.division||0,'Cumulative divisions do not depend on which chunks were viewed')}}
function renderedBodies(){return JSON.parse(vm.runInContext('JSON.stringify(points.map(p=>({id:p.c.cell_id,x:p.x,y:p.y,r:p.r,polygon:p.polygon,worldRadius:p.worldRadius,measured:p.measured,synthetic:p.synthetic})))',realm))}
function width(p){return Math.max(...p.polygon.map(v=>v[0]))-Math.min(...p.polygon.map(v=>v[0]))}
function close(actual,expected,label){assert(Math.abs(actual-expected)<1e-7,label+`: ${actual} versus ${expected}`)}
function bodyFill(p){return mapFills.find(fill=>fill.vertices.length===p.polygon.length&&fill.vertices[0]?.every((v,i)=>Math.abs(v-p.polygon[0][i])<1e-7))}
function nucleus(p){return mapFills.find(fill=>fill.color==='#0a0e1a'&&fill.arcs.some(a=>Math.abs(a.x-p.x)<1e-7&&Math.abs(a.y-p.y)<1e-7))?.arcs[0]}
function checkOriginalGeometry(){
 const original=fs.readFileSync(path.join(__dirname,'../src/virtualtissue_xenium/engine/web/player_template.html'),'utf8');
 const start=original.indexOf('function powerCells(cells){'),end=original.indexOf('// ---------------------------------------------------------------- drawing',start);
 assert(start>=0&&end>start,'The authoritative original cell renderer is available');
 const oldRealm=vm.createContext({});vm.runInContext(original.slice(start,end),oldRealm);
 for(const sites of [[{x:0,y:0,r:8}],
                     [{x:0,y:0,r:8},{x:10,y:0,r:8}],
                     [{x:0,y:0,r:12},{x:11,y:3,r:5},{x:-8,y:6,r:7},{x:0,y:0,r:3}]]){
  oldRealm.sites=sites;realm.sites=sites;
  assert.equal(vm.runInContext('JSON.stringify(powerCells(sites))',realm),
               vm.runInContext('JSON.stringify(powerCells(sites))',oldRealm),
               'Isolated, touching and unequal-sized silhouettes exactly match the original renderer');
 }
}
async function checkShapes(){
 await waitFor(()=>ids.canvas.dataset.frameStep==='0');
 checkOriginalGeometry();
 assert.equal(ids.shapeMode.value,'cells','Cell silhouettes are the initial view');
 let rendered=Object.fromEntries(renderedBodies().map(p=>[p.id,p]));
 const solo=rendered.cell,cd8=rendered.cd8,a=rendered.pair_a,b=rendered.pair_b,provided=rendered.provided;
 assert.equal(solo.polygon.length,40,'An isolated cell uses the original40-sided silhouette');
 assert(a.polygon.length>=3&&a.polygon.length<40,'Touching cell outlines are clipped to their shared boundary');
 const midpoint=(a.x+b.x)/2;
 assert(Math.max(...a.polygon.map(p=>p[0]))<midpoint);
 assert(Math.min(...b.polygon.map(p=>p[0]))>midpoint,'Clipped neighboring bodies do not overlap');
 close(solo.worldRadius,8,'Saved malignant radius is preserved');close(cd8.worldRadius,5,'Saved Tcell radius is preserved');
 close(solo.r/cd8.r,8/5,'Different lineage sizes share the same world-to-screen scale');
 close(rendered.default_radius.worldRadius,5,'A cell without a saved radius uses the ontology fallback');
 assert.equal(bodyFill(solo).color,'#ff8a5c','Tumor silhouette uses the old coral palette');
 assert.equal(bodyFill(cd8).color,'#38d6cf','CD8 silhouette uses the old turquoise palette');
 assert(nucleus(solo)&&nucleus(cd8),'Visible cell bodies have dark nuclei');
 assert(nucleus(solo).r/solo.r>.2&&nucleus(solo).r/solo.r<.35,'Nuclei scale with cell bodies');
 assert.equal(provided.measured,true);assert.equal(provided.synthetic,false);
 assert.equal(provided.polygon.length,5,'Supplied square is retained instead of becoming a40-sided synthetic body');
 const measured=data.replay?JSON.parse(vm.runInContext('JSON.stringify(currentFrame.cells.provided)',realm)):data.frames[0].cells.provided;
 const pixelScale=solo.r/8;
 const checkProvided=(p,scale)=>{for(let i=0;i<p.polygon.length;i++){
  const raw=measured.t0_geometry.polygon_um[i],anchor=measured.t0_geometry.anchor_um;
  close(p.polygon[i][0]-p.x,(raw[0]-anchor[0])*scale,'Provided x vertex is exact');
  close(p.polygon[i][1]-p.y,(raw[1]-anchor[1])*scale,'Provided y vertex is exact');
 }};
 checkProvided(provided,pixelScale);
 const before={radius:solo.r,width:width(solo),nucleus:nucleus(solo).r};
 ids.canvas.onwheel({preventDefault(){},clientX:300,clientY:280,deltaY:-Math.log(2)/.0015});
 rendered=Object.fromEntries(renderedBodies().map(p=>[p.id,p]));
 close(rendered.cell.r/before.radius,2,'Zoom doubles the cell radius');
 close(width(rendered.cell)/before.width,2,'Zoom doubles the full silhouette');
 close(nucleus(rendered.cell).r/before.nucleus,2,'Zoom doubles the nucleus');
 ids.shapeMode.value='centroids';ids.shapeMode.onchange();
 assert(renderedBodies().every(p=>p.polygon===null),'The explicit centroid view retains simple markers');
 assert.equal(ids.canvas.dataset.providedOutlines,'0');
 ids.shapeMode.value='provided';ids.shapeMode.onchange();
 rendered=Object.fromEntries(renderedBodies().map(p=>[p.id,p]));
 assert.equal(rendered.cell.polygon,null);assert.equal(rendered.provided.measured,true);
 ids.shapeMode.value='cells';ids.shapeMode.onchange();
 assert.equal(ids.canvas.dataset.illustrativeBodies,'5');assert.equal(ids.canvas.dataset.providedOutlines,'1');
 ids.slider.value='1';await ids.slider.oninput();
 rendered=Object.fromEntries(renderedBodies().map(p=>[p.id,p]));
 const midpointCell=rendered.provided,scale=rendered.cell.r/8;
 const sourceFill=mapFills.find(fill=>fill.vertices.length===5&&Math.abs(fill.vertices[0][0]-(midpointCell.x-4*scale))<1e-7);
 const syntheticFill=mapFills.find(fill=>fill.vertices.length===40&&Math.abs(fill.vertices.reduce((sum,p)=>sum+p[0],0)/40-midpointCell.x)<1e-7);
 assert(sourceFill&&syntheticFill,'At10minutes both measured and original silhouettes are visible');
 close(sourceFill.alpha,syntheticFill.alpha,'The original20minute transition reaches its midpoint');
 assert(sourceFill.alpha>0&&sourceFill.alpha<bodyFill(rendered.cell).alpha);
 for(const step of [2,3]){
  ids.slider.value=String(step);await ids.slider.oninput();
  rendered=Object.fromEntries(renderedBodies().map(p=>[p.id,p]));
  assert.equal(rendered.provided.measured,false);assert.equal(rendered.provided.synthetic,true);
  assert.equal(rendered.provided.polygon.length,40,'After20minutes original playback uses the old cell bodies');
 }
 ids.shapeMode.value='provided';ids.shapeMode.onchange();
 rendered=Object.fromEntries(renderedBodies().map(p=>[p.id,p]));
 assert.equal(rendered.provided.measured,true);assert.equal(rendered.provided.polygon.length,5);
 checkProvided(rendered.provided,scale);
 assert.equal(vm.runInContext('JSON.stringify(currentFrame.cells.provided.t0_geometry)',realm),JSON.stringify(measured.t0_geometry),'The visual transition never changes the source geometry');
 console.log(JSON.stringify({passed:true,checked:['original-renderer-equivalence','organic-default','40-sided-isolated-cell','neighbor-clipping','saved-world-radii','old-palette','visible-nuclei','exact-provided-outlines','zoom','centroid-mode','original20minute-transition']}));
}
(async()=>{
 if(process.argv[3]==='shapes'){await checkShapes();return}
 await waitFor(()=>ids.canvas.dataset.frameStep==='0');assert.equal(requested.length,1);assert.equal(body.dataset.replayCachedChunks,'1');
 if(process.argv[3]==='provided'){ids.shapeMode.value='provided';ids.shapeMode.onchange()}
 assert(chartSeriesLengths.includes(data.frames.length),'The population chart includes the complete saved history while only one chunk is loaded');
 await seek(55);await seek(110);assert.equal(body.dataset.replayCachedChunks,'2');assert.equal(requested.length,3);
 await seek(0);assert.equal(requested.length,4,'Evicted first chunk must load again');
 const before=requested.length;const jobs=[];for(const step of [51,101,150]){ids.slider.value=String(step);jobs.push(ids.slider.oninput())}await Promise.all(jobs);
 assert.equal(ids.canvas.dataset.frameStep,'150');assert.equal(requested.length,before+1,'Obsolete pending seeks must not load intermediate chunks');
 await seek(100);ids.canvas.onclick({clientX:points[0].x,clientY:points[0].y});assert.equal(JSON.parse(ids.detail.textContent).cell_id,'cell');
 assert.equal(ids.tracks.checked,false);const previousStrokes=strokes;ids.tracks.checked=true;ids.tracks.onchange();assert(strokes>=previousStrokes+2,'Tracks add a movement segment in addition to the selected-cell outline');
 ids.play.onclick();await waitFor(()=>ids.canvas.dataset.frameStep==='101');ids.play.onclick();assert.equal(ids.play.textContent,'Play');
 console.log(JSON.stringify({passed:true,loadedRequests:requested.length,cacheChunks:Number(body.dataset.replayCachedChunks),frame:ids.canvas.dataset.frameStep,polygonDraws,polygonVertexCounts,checked:['initial-load','seek','LRU-eviction','rapid-seek','selection','recent-tracks','play-pause']}));
})().catch(error=>{console.error(error);process.exitCode=1});
