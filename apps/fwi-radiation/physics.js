/* Isotropic Born kernels: homogeneous plane waves, complex scalar slowness, e^-iwt. */
const C=(r=0,i=0)=>[r,i], add=(a,b)=>[a[0]+b[0],a[1]+b[1]], mul=(a,b)=>[a[0]*b[0]-a[1]*b[1],a[0]*b[1]+a[1]*b[0]], scale=(a,x)=>[a[0]*x,a[1]*x], dot=(a,b)=>a.reduce((s,x,i)=>s+x*b[i],0);
const div=(a,b)=>{let d=b[0]*b[0]+b[1]*b[1];return[(a[0]*b[0]+a[1]*b[1])/d,(a[1]*b[0]-a[0]*b[1])/d];};
const sqrtc=a=>{let h=Math.hypot(...a);return[Math.sqrt(Math.max(0,(h+a[0])/2)),(a[1]<0?-1:1)*Math.sqrt(Math.max(0,(h-a[0])/2))];};
function moduli(s){let visco=s.visco||s.basis.endsWith('-q'),rho=s.rho,M=rho*s.vp*s.vp,U=s.elastic?rho*s.vs*s.vs:0,L=C(visco?2/Math.PI*Math.log(s.f/s.fref):0,visco?-1:0);return{M,U,L,mc:add(C(M),scale(L,M/s.qp)),uc:add(C(U),scale(L,U/s.qs))};}
function lowLossSize(s){if(!(s.visco||s.basis.endsWith('-q')))return 0;let l=Math.abs(2/Math.PI*Math.log(s.f/s.fref));return Math.max(l/s.qp,s.elastic?l/s.qs:0);}
function validateMedium(s){if(!(s.vp>0&&s.rho>0&&(!s.elastic||s.vs>0)&&s.f>0&&s.fref>0))return 'Use positive velocities, density and frequencies.';let {M,U,L,mc,uc}=moduli(s);if(M<=4*U/3)return 'Choose Vp² > 4Vs²/3 for a positive bulk modulus.';if(s.visco||s.basis.endsWith('-q')){if(!(s.qp>=20&&(!s.elastic||s.qs>=20)))return 'Use Q ≥ 20 for the low-loss model.';if(s.elastic&&M/s.qp+1e-12<4*U/(3*s.qs))return `Non-passive bulk attenuation: Qp must be ≤ ${(3/4*(s.vp/s.vs)**2*s.qs).toFixed(2)} for these Vp, Vs and Qs values. Reduce Qp or increase Qs.`;if(mc[0]<=4*uc[0]/3||uc[0]<0)return 'The real bulk or shear modulus is non-positive at this frequency.';}return '';}
function gramWeight(s,mode,angle){let {M,mc,uc}=moduli(s),out=mode[1]==='P'?mc:uc,w=s.modeWeight==='displacement'?M*M/(out[0]**2+out[1]**2):1;return w*(s.angleMeasure==='solid'?Math.max(0,Math.sin(angle*Math.PI/180)):1);}
function parameters(s){let b=s.basis,a;
if(b==='velocity'||b==='velocity-q')a=['ln Vp',...(s.elastic?['ln Vs']:[]),...(b==='velocity'?['ln ρ']:[])];
else if(b==='impedance'||b==='impedance-q')a=['ln Ip',...(s.elastic?['ln Is']:[]),...(b==='impedance'?['ln ρ']:[])];
else if(b==='velocity-impedance')a=['ln Vp',...(s.elastic?['ln Vs']:[]),'ln Ip'];
else if(b==='velocity-s-impedance')a=['ln Vp','ln Vs','ln Is'];
else if(b==='modulus-velocity')a=['ln K','ln Vp'];
else if(b==='modulus-impedance')a=['ln K','ln Ip'];
else a=s.elastic?['ln K','ln μ','ln ρ']:['ln K','ln ρ'];
if(s.visco||b.endsWith('-q'))a.push(s.qcoord==='quality'?'ln Qp':'qP = 1/Qp',...(s.elastic?[s.qcoord==='quality'?'ln Qs':'qS = 1/Qs']:[]));return a;}
function kernel(s,angle,mode,p){let vp=s.vp,vs=s.elastic?s.vs:0,rho=s.rho,M=rho*vp*vp,U=rho*vs*vs,K=M-4*U/3;let L=C((s.visco||s.basis.endsWith('-q'))?2/Math.PI*Math.log(s.f/s.fref):0,(s.visco||s.basis.endsWith('-q'))?-1:0),mc=add(C(M),scale(L,M/s.qp)),uc=add(C(U),scale(L,U/s.qs));let dm=C(),du=C(),dr=0;
if(p==='qP = 1/Qp')dm=scale(L,M*.01);
else if(p==='qS = 1/Qs')du=scale(L,U*.01);
else if(p==='ln Qp')dm=scale(L,-M/s.qp);
else if(p==='ln Qs')du=scale(L,-U/s.qs);
else if(s.basis==='modulus'){if(p==='ln K')dm=scale(mc,K/M);if(p==='ln μ'){du=uc;dm=scale(mc,4*U/(3*M));}if(p==='ln ρ')dr=rho;}
else if(s.basis==='modulus-velocity'){if(p==='ln K'){dm=mc;dr=rho;}if(p==='ln Vp')dr=-2*rho;}
else if(s.basis==='modulus-impedance'){if(p==='ln K'){dm=mc;dr=-rho;}if(p==='ln Ip')dr=2*rho;}
else {let v=0,w=0,r=0,b=s.basis;
if(b==='velocity'||b==='velocity-q'){v=+(p==='ln Vp');w=+(p==='ln Vs');r=+(p==='ln ρ');}
if(b==='impedance'||b==='impedance-q'){r=+(p==='ln ρ');v=+(p==='ln Ip')-r;w=+(p==='ln Is')-r;}
if(b==='velocity-impedance'){v=+(p==='ln Vp');w=+(p==='ln Vs');r=+(p==='ln Ip')-v;}
if(b==='velocity-s-impedance'){v=+(p==='ln Vp');w=+(p==='ln Vs');r=+(p==='ln Is')-w;}
dm=scale(mc,r+2*v);du=scale(uc,r+2*w);dr=rho*r;}

let t=angle*Math.PI/180,ni=[0,0,1],no=[Math.sin(t),0,Math.cos(t)],ei=mode[0]==='P'?ni:mode==='HH'?[0,1,0]:[1,0,0],eo=mode[1]==='P'?no:mode==='HH'?[0,1,0]:[Math.cos(t),0,-Math.sin(t)],vi=mode[0]==='P'?vp:vs,vo=mode[1]==='P'?vp:vs;
let dl=add(dm,scale(du,-2)),a=dot(eo,no)*dot(ni,ei),b=dot(eo,ei)*dot(no,ni)+dot(eo,ni)*dot(no,ei),Mi=mode[0]==='P'?mc:uc,Mo=mode[1]==='P'?mc:uc,kk=div(C(rho),mul(sqrtc(Mi),sqrtc(Mo)));
return scale(add(C(dr*dot(eo,ei)),scale(mul(add(scale(dl,a),scale(du,b)),kk),-1)),1/rho);}

function compute(s){let error=validateMedium(s);if(error)throw new Error(error);let ps=parameters(s),modes=s.elastic?s.modes:['PP'],plotMode=modes.includes(s.plotMode)?s.plotMode:modes[0],usedModes=s.matrixScope==='stack'?modes:[plotMode],rows=[];for(let a=0;a<=180;a++)for(let m of modes)rows.push({a,m,k:ps.map(p=>kernel(s,a,m,p)),weight:gramWeight(s,m,a),used:a>=s.amin&&a<=s.amax&&usedModes.includes(m)});let selected=rows.filter(r=>r.used),h=ps.map((p,i)=>ps.map((q,j)=>selected.reduce((sum,r)=>sum+r.weight*(s.metric==='complex'?dot(r.k[i],r.k[j]):Math.hypot(...r.k[i])*Math.hypot(...r.k[j])),0)));let norms=ps.map((p,i)=>Math.sqrt(Math.max(0,h[i][i]))),corr=h.map((r,i)=>r.map((v,j)=>norms[i]*norms[j]>1e-12?Math.max(-1,Math.min(1,v/(norms[i]*norms[j]))):null));return{ps,rows,h,norms,corr,usedModes,plotMode};}
if(typeof module!=='undefined')module.exports={parameters,kernel,compute,moduli,validateMedium,gramWeight,lowLossSize};
