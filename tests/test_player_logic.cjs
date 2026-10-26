// Offline JavaScript regression checks. This is not a browser/rendering test.
const fs = require('node:fs');
const path = require('node:path');
const vm = require('node:vm');
const assert = require('node:assert/strict');
const html = fs.readFileSync(path.join(__dirname, '../src/virtualtissue_xenium/engine/web/player_template.html'), 'utf8');
const source = html.slice(html.indexOf('"use strict";'), html.lastIndexOf('load().then('));
const names = {oxygen:'Oxygen', ifng:'IFN-γ-like',cxcl9:'CXCL9/10-like',tgfb:'TGF-β-like',
  glucose:'Glucose',lactate:'Lactate',vegf:'VEGF-like',matrix:'Matrix',ne:'NE-like',cgrp:'CGRP-like'};
const zeros = Buffer.alloc(9).toString('base64');
function payload() {
  const fields = t => ({t,stats:Object.fromEntries(Object.keys(names).map(k=>[k,{min:0,max:0}])),
    ...Object.fromEntries(Object.keys(names).map(k=>[k,zeros]))});
  return {patch:'example',sample_name:'Example sample',provider:'local_score_rules',scenario:'declared_recognition',
    seed:21,stop:'horizon_reached',paid_calls:0,scope:'Test fixture',signatures:[],
    actions:['','CYCLING'],events_kinds:['division','apoptosis_committed','intact_death','clearance','corpse_degraded','hypoxic_necrosis'],
    cells:{1:{l:'Tumor',o:'roi',b:0,u:null},2:{l:'T_CD8_cytotoxic',o:'roi',b:0,u:null}},
    frames:[{t:0,ids:[1,2],x:[0,10],y:[0,10],r:[2,2],a:[0,0],tg:[0,0],v:[0,0],cp:[-1,-1],cb:[0,0]},
            {t:60,ids:[2],x:[10],y:[10],r:[2],a:[0],tg:[0],v:[0],cp:[-1],cb:[0]}],
    events:[[60,5,1,0,0,0,false]],decisions:{},fields:[fields(0),fields(60)],auto_field:'oxygen',
    field_grid:{origin:[0,0],spacing:5,shape:[3,3],names,scales:Object.fromEntries(Object.keys(names).map(k=>[k,1])),statuses:{ifng:'disabled'}},
    population_limit:{maximum:2,initial_headroom:0,reached_min:0}};
}

for (const [width,height] of [[1400,900],[1920,1080]]) {
  const text=[],elements=new Map();
  const ctx=new Proxy({measureText:s=>({width:s.length*7}),fillText:(...args)=>text.push(args),
    createImageData:(w,h)=>({data:new Uint8ClampedArray(w*h*4)}),
    createLinearGradient:()=>({addColorStop(){}}),createRadialGradient:()=>({addColorStop(){}})}, {get:(o,k)=>k in o?o[k]:(()=>{})});
  function element(id) {
    if(!elements.has(id))elements.set(id,{value:id==='field'?'auto':id==='colorMode'?'type':'',checked:true,style:{},
      textContent:'',innerHTML:'',getContext:()=>ctx,querySelectorAll:()=>[],setAttribute(){},addEventListener(){},
      getBoundingClientRect:()=>({width:width-360,height:height-220,top:0,bottom:height-220})});
    return elements.get(id);
  }
  const sandbox=vm.createContext({document:{getElementById:element,createElement:()=>element('temporaryCanvas'),
      querySelector:()=>({getBoundingClientRect:()=>({bottom:232})}),querySelectorAll:()=>[],addEventListener(){}},
    window:{devicePixelRatio:1,addEventListener(){}},performance:{now:()=>0},
    requestAnimationFrame(){},Path2D:class{moveTo(){}lineTo(){}closePath(){}arc(){}},
    atob:s=>Buffer.from(s,'base64').toString('binary'), fixture:payload()});
  vm.runInContext(source,sandbox);
  vm.runInContext('P=fixture;init();T=60;draw();updatePanel();',sandbox);
  assert.equal(element('tDeath').textContent,1);
  assert.match(element('fieldSummary').textContent,/Oxygen = 0.*anoxia/);
  assert.match(element('highlights').innerHTML,/Death from hypoxia/);
  assert.match(element('ceilingStatus').textContent,/computational limit/);
  assert.ok(text.some(x=>String(x[0]).includes('Death from hypoxia')));
  const before=text.length;
  vm.runInContext('view.x=0;view.y=500;view.s=1;drawHighlights();',sandbox);
  assert.ok(text.slice(before).some(x=>String(x[0]).includes('(map edge)')));
  vm.runInContext('T=0;draw();updatePanel();',sandbox);
  assert.equal(element('tDeath').textContent,0);
  vm.runInContext('T=60;draw();updatePanel();',sandbox);
  assert.equal(element('tDeath').textContent,1);
  element('field').value='ifng';vm.runInContext('draw();',sandbox);
  assert.match(element('fieldSummary').textContent,/Not produced/);
  assert.doesNotMatch(element('fieldStatus').innerHTML,/scale-ends/);
  element('field').value='cxcl9';vm.runInContext('P.field_grid.scale_has_signal={cxcl9:false};draw();',sandbox);
  assert.doesNotMatch(element('fieldStatus').innerHTML,/scale-ends/);
  element('field').value='ne';vm.runInContext('delete P.fields[1].ne;fieldFrameCache=null;draw();',sandbox);
  assert.match(element('fieldSummary').textContent,/not saved/);
  element('field').value='auto';vm.runInContext('draw();',sandbox);
  assert.match(element('fieldSummary').textContent,/Oxygen/);
  console.log(`PASS: ${width}x${height} drawing logic, zero oxygen, removed-cell event, backward seek, scenario and missing-field states`);
}
