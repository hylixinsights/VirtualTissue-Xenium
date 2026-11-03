// A DOM test double executes the actual lazy player and its generated sidecars.
// This tests loader/cache/seek behavior, not browser visual rendering.
const assert=require('node:assert/strict'),fs=require('node:fs'),path=require('node:path'),vm=require('node:vm');
const folder=process.argv[2],html=fs.readFileSync(path.join(folder,'replay.html'),'utf8');
const encoded=html.match(/<script id="run-data" type="application\/json">([\s\S]*?)<\/script>/)[1];
const player=html.split('</script><script>')[1].split('</script>')[0];
let points=[],strokes=0,polygonDraws=0,currentPath=[];const requested=[];
const context2d={setTransform(){},clearRect(){points=[]},beginPath(){currentPath=[]},arc(x,y,r){points.push({x,y,r})},fill(){},stroke(){strokes++},moveTo(x,y){currentPath=[[x,y]]},lineTo(x,y){currentPath.push([x,y])},closePath(){polygonDraws++;points.push({x:currentPath.reduce((s,p)=>s+p[0],0)/currentPath.length,y:currentPath.reduce((s,p)=>s+p[1],0)/currentPath.length,r:5})},setLineDash(){}};
function element(id=''){return {id,dataset:{},style:{},textContent:'',value:'0',checked:false,disabled:false,append(){},remove(){},getBoundingClientRect(){return {width:600,height:560,left:0,top:0}},getContext(){return context2d}}}
const ids=Object.fromEntries(['run-data','slider','canvas','play','tracks','warning','summary','legend','time','stats','detail'].map(id=>[id,element(id)]));
ids['run-data'].textContent=encoded;const body=element('body'),heading=element('h1');
let realm;
const document={
 body,title:'',getElementById:id=>ids[id],querySelector:selector=>heading,
 createElement:tag=>element(tag),createTextNode:text=>({textContent:text}),
 head:{append(script){
  requested.push(script.src);
  setImmediate(()=>{
   try{vm.runInContext(fs.readFileSync(path.join(folder,script.src),'utf8'),realm);if(script.onload)script.onload()}
   catch(error){if(script.onerror)script.onerror(error)}
  });
 }}
};
realm=vm.createContext({document,devicePixelRatio:1,setTimeout,clearTimeout,console});
realm.window=realm;realm.window.addEventListener=()=>{};
vm.runInContext(player,realm);
async function waitFor(predicate){for(let i=0;i<100;i++){if(predicate())return;await new Promise(resolve=>setTimeout(resolve,5))}throw new Error('Timed out waiting for replay state: '+JSON.stringify(body.dataset))}
async function seek(step){ids.slider.value=String(step);await ids.slider.oninput();assert.equal(ids.canvas.dataset.frameStep,String(step));assert.equal(body.dataset.replayStatus,'ready');assert(Number(body.dataset.replayCachedChunks)<=2);assert(Number(body.dataset.replayCachedFrames)<=100)}
(async()=>{
 await waitFor(()=>ids.canvas.dataset.frameStep==='0');assert.equal(requested.length,1);assert.equal(body.dataset.replayCachedChunks,'1');
 await seek(55);await seek(110);assert.equal(body.dataset.replayCachedChunks,'2');assert.equal(requested.length,3);
 await seek(0);assert.equal(requested.length,4,'Evicted first chunk must load again');
 const before=requested.length;const jobs=[];for(const step of [51,101,150]){ids.slider.value=String(step);jobs.push(ids.slider.oninput())}await Promise.all(jobs);
 assert.equal(ids.canvas.dataset.frameStep,'150');assert.equal(requested.length,before+1,'Obsolete pending seeks must not load intermediate chunks');
 await seek(100);ids.canvas.onclick({clientX:points[0].x,clientY:points[0].y});assert.equal(JSON.parse(ids.detail.textContent).cell_id,'cell');
 assert.equal(ids.tracks.checked,false);const previousStrokes=strokes;ids.tracks.checked=true;ids.tracks.onchange();assert(strokes>=previousStrokes+2,'Tracks add a movement segment in addition to the selected-cell outline');
 ids.play.onclick();await waitFor(()=>ids.canvas.dataset.frameStep==='101');ids.play.onclick();assert.equal(ids.play.textContent,'Play');
 console.log(JSON.stringify({passed:true,loadedRequests:requested.length,cacheChunks:Number(body.dataset.replayCachedChunks),frame:ids.canvas.dataset.frameStep,polygonDraws,checked:['initial-load','seek','LRU-eviction','rapid-seek','selection','recent-tracks','play-pause']}));
})().catch(error=>{console.error(error);process.exitCode=1});
