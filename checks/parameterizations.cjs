const assert=require('node:assert/strict');const{kernel,parameters}=require('../apps/fwi-radiation/physics.js');
const dot=(a,b)=>a.reduce((v,x,i)=>v+x*b[i],0);let count=0;
for(let elastic of [false,true])for(let visco of [false,true])for(let qcoord of ['quality','inverse'])for(let basis of ['velocity','impedance','velocity-impedance','modulus','velocity-q','impedance-q',...(elastic?['velocity-s-impedance']:['modulus-velocity','modulus-impedance'])]){
let s={elastic,visco:visco||basis.endsWith('-q'),qcoord,basis,vp:3.2,vs:1.4,rho:2.2,qp:80,qs:55,f:17,fref:10};
function raw(p,e){let v={Vp:s.vp,Vs:elastic?s.vs:0,'ρ':s.rho,Ip:s.rho*s.vp,Is:s.rho*s.vs,K:s.rho*(s.vp**2-4*(elastic?s.vs**2:0)/3),'μ':elastic?s.rho*s.vs**2:0,Qp:s.qp,Qs:s.qs,qP:1/s.qp,qS:1/s.qs};if(p.startsWith('ln '))v[p.slice(3)]*=Math.exp(e);else v[p[1]==='P'?'qP':'qS']+=.01*e;
let rho=v['ρ'],vp=v.Vp,vs=elastic?v.Vs:0,M,U;
if(basis==='velocity-impedance')rho=v.Ip/vp;
if(basis==='velocity-s-impedance')rho=v.Is/vs;
if(basis==='impedance'||basis==='impedance-q'){vp=v.Ip/rho;if(elastic)vs=v.Is/rho;}
if(basis==='modulus'){M=v.K+4*v['μ']/3;U=v['μ'];}
else if(basis==='modulus-velocity'){rho=v.K/vp**2;M=v.K;U=0;}
else if(basis==='modulus-impedance'){rho=v.Ip**2/v.K;M=v.K;U=0;}
else {M=rho*vp**2;U=rho*vs**2;}
let qp=qcoord==='quality'?1/v.Qp:v.qP,qs=qcoord==='quality'?1/v.Qs:v.qS,L=2/Math.PI*Math.log(s.f/s.fref);return{rho,M:[M*(1+(s.visco?qp*L:0)),s.visco?-M*qp:0],U:[U*(1+(s.visco?qs*L:0)),s.visco?-U*qs:0]};}
for(let p of parameters(s))for(let a of [0,37,90,151,180])for(let mode of elastic?['PP','PS','SP','SS','HH']:['PP']){let t=a*Math.PI/180,ni=[0,0,1],no=[Math.sin(t),0,Math.cos(t)],ei=mode[0]==='P'?ni:mode==='HH'?[0,1,0]:[1,0,0],eo=mode[1]==='P'?no:mode==='HH'?[0,1,0]:[Math.cos(t),0,-Math.sin(t)],A=dot(eo,no)*dot(ni,ei),D=dot(eo,ei)*dot(no,ni)+dot(eo,ni)*dot(no,ei),vv=(mode[0]==='P'?s.vp:s.vs)*(mode[1]==='P'?s.vp:s.vs),h=1e-5,plus=raw(p,h),minus=raw(p,-h),k=kernel(s,a,mode,p);
let base=raw(p,0),mi=mode[0]==='P'?base.M:base.U,mo=mode[1]==='P'?base.M:base.U;
// Independently construct slownesses with polar-form square roots.
function slow(z){let mag=Math.hypot(...z),phi=Math.atan2(z[1],z[0]);return[Math.sqrt(s.rho/mag)*Math.cos(-phi/2),Math.sqrt(s.rho/mag)*Math.sin(-phi/2)];}
let si=slow(mi),so=slow(mo),kk=[si[0]*so[0]-si[1]*so[1],si[0]*so[1]+si[1]*so[0]],delta=[0,1].map(j=>(((plus.M[j]-minus.M[j])-2*(plus.U[j]-minus.U[j]))*A+(plus.U[j]-minus.U[j])*D)/(2*h)),stiff=[delta[0]*kk[0]-delta[1]*kk[1],delta[0]*kk[1]+delta[1]*kk[0]];
for(let j=0;j<2;j++){let dr=j===0?(plus.rho-minus.rho)/(2*h):0,fd=(dr*dot(eo,ei)-stiff[j])/s.rho;assert.ok(Math.abs(fd-k[j])<1e-7,JSON.stringify({basis,elastic,visco,qcoord,p,a,mode,fd,k:k[j]}));}count++;}}
console.log(`PASS: ${count} centered finite-difference checks across all supported bases, media, Q coordinates and modes.`);
