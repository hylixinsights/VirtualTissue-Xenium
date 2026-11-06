// A DOM test double executes the actual lazy player and its generated sidecars.
// This tests loader/cache/seek behavior, not browser visual rendering.
const assert=require('node:assert/strict'),fs=require('node:fs'),path=require('node:path'),vm=require('node:vm');
const folder=process.argv[2],html=fs.readFileSync(path.join(folder,'replay.html'),'utf8');
const encoded=html.match(/<script id="run-data" type="application\/json">([\s\S]*?)<\/script>/)[1];
const data=JSON.parse(encoded);
const player=html.split('</script><script>')[1].split('</script>')[0];
let points=[],strokes=0,polygonDraws=0,currentPath=[];const requested=[],chartSeriesLengths=[];
function context2d(id){
 let drawing=[];
 return {
  setTransform(){},clearRect(){if(id==='canvas')points=[]},
  beginPath(){drawing=[];if(id==='canvas')currentPath=[]},
  arc(x,y,r){if(id==='canvas')points.push({x,y,r})},
  fill(){},fillRect(){},fillText(){},
  stroke(){if(id==='canvas')strokes++;if(id==='chart')chartSeriesLengths.push(drawing.length)},
  moveTo(x,y){drawing=[[x,y]];if(id==='canvas')currentPath=[[x,y]]},
  lineTo(x,y){drawing.push([x,y]);if(id==='canvas')currentPath.push([x,y])},
  closePath(){if(id==='canvas'&&currentPath.length){polygonDraws++;points.push({x:currentPath.reduce((s,p)=>s+p[0],0)/currentPath.length,y:currentPath.reduce((s,p)=>s+p[1],0)/currentPath.length,r:5})}},
  setLineDash(){},createRadialGradient(){return {addColorStop(){}}},createLinearGradient(){return {addColorStop(){}}},
  measureText(t){return {width:String(t).length*6}},roundRect(){},save(){},restore(){},translate(){},scale(){}
 };
}
function element(id=''){const ctx=context2d(id);return {id,dataset:{},style:{},textContent:'',innerHTML:'',value:'0',checked:false,disabled:false,children:[],append(...items){this.children.push(...items)},appendChild(item){this.append(item)},remove(){},setAttribute(k,v){this[k]=v},getAttribute(k){return this[k]},addEventListener(k,fn){this['on'+k]=fn},classList:{add(){},remove(){},toggle(){}},querySelectorAll(){return this.children.filter(x=>x.tagName==='BUTTON')},getBoundingClientRect(){return {width:600,height:id==='chart'?110:560,left:0,top:0,bottom:560}},getContext(){return ctx}}}
const ids={};
for(const tag of html.matchAll(/<[^>]+\bid="([^"]+)"[^>]*>/g)){
 const e=element(tag[1]);e.checked=/\bchecked\b/.test(tag[0]);const value=tag[0].match(/\bvalue="([^"]*)"/);if(value)e.value=value[1];ids[e.id]=e;
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
(async()=>{
 await waitFor(()=>ids.canvas.dataset.frameStep==='0');assert.equal(requested.length,1);assert.equal(body.dataset.replayCachedChunks,'1');
 assert(chartSeriesLengths.includes(data.frames.length),'The population chart includes the complete saved history while only one chunk is loaded');
 await seek(55);await seek(110);assert.equal(body.dataset.replayCachedChunks,'2');assert.equal(requested.length,3);
 await seek(0);assert.equal(requested.length,4,'Evicted first chunk must load again');
 const before=requested.length;const jobs=[];for(const step of [51,101,150]){ids.slider.value=String(step);jobs.push(ids.slider.oninput())}await Promise.all(jobs);
 assert.equal(ids.canvas.dataset.frameStep,'150');assert.equal(requested.length,before+1,'Obsolete pending seeks must not load intermediate chunks');
 await seek(100);ids.canvas.onclick({clientX:points[0].x,clientY:points[0].y});assert.equal(JSON.parse(ids.detail.textContent).cell_id,'cell');
 assert.equal(ids.tracks.checked,false);const previousStrokes=strokes;ids.tracks.checked=true;ids.tracks.onchange();assert(strokes>=previousStrokes+2,'Tracks add a movement segment in addition to the selected-cell outline');
 ids.play.onclick();await waitFor(()=>ids.canvas.dataset.frameStep==='101');ids.play.onclick();assert.equal(ids.play.textContent,'Play');
 console.log(JSON.stringify({passed:true,loadedRequests:requested.length,cacheChunks:Number(body.dataset.replayCachedChunks),frame:ids.canvas.dataset.frameStep,polygonDraws,checked:['initial-load','seek','LRU-eviction','rapid-seek','selection','recent-tracks','play-pause']}));
})().catch(error=>{console.error(error);process.exitCode=1});
