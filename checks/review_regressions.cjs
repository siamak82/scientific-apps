const assert=require('node:assert/strict');const{compute,kernel,validateMedium,gramWeight}=require('../apps/fwi-radiation/physics.js');
const base={elastic:true,visco:false,basis:'velocity',qcoord:'quality',vp:3,vs:1.5,rho:2,qp:80,qs:50,f:10,fref:10,amin:0,amax:180,metric:'complex',modes:['PP','PS'],matrixScope:'stack',modeWeight:'equal',angleMeasure:'section',plotMode:'PP'};
const uniform=compute(base),displacement=compute({...base,modeWeight:'displacement'});
assert.ok(uniform.corr[1][2]<0&&displacement.corr[1][2]>0);
assert.equal(gramWeight({...base,modeWeight:'displacement'},'PS',90),16);
const one=compute({...base,matrixScope:'single'}),pp=compute({...base,modes:['PP']});assert.deepEqual(one.corr,pp.corr);
const ps=compute({...base,matrixScope:'single',plotMode:'PS'});assert.deepEqual(ps.usedModes,['PS']);assert.equal(ps.norms[0],0);
const bad={...base,visco:true,qp:200};assert.match(validateMedium(bad),/Non-passive/);assert.throws(()=>compute(bad),/Non-passive/);assert.equal(validateMedium({...bad,qp:150}),'');
const ac={...base,elastic:false,visco:true,modes:['PP'],qp:20};for(let a of [0,45,90,180]){assert.ok(Math.abs(kernel(ac,a,'PP','ln Vp')[1])<1e-12);assert.ok(Math.abs(kernel(ac,a,'PP','ln ρ')[1])<1e-12);}
let section=compute(base),solid=compute({...base,angleMeasure:'solid'});assert.ok(Math.abs(section.corr[0][1]-solid.corr[0][1])>.01);assert.equal(gramWeight({...base,angleMeasure:'solid'},'PP',0),0);
console.log('Uniform PP+PS correlations:',[uniform.corr[0][1],uniform.corr[0][2],uniform.corr[1][2]].map(x=>x.toFixed(3)).join(', '));
console.log('Displacement PP+PS correlations:',[displacement.corr[0][1],displacement.corr[0][2],displacement.corr[1][2]].map(x=>x.toFixed(3)).join(', '));
console.log('PASS: mode factors, single-mode default, measure selection, passivity rejection and complex-kinematics limits.');
