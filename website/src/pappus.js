// pappus.js - compact prototype
// Basic helpers
const canvas = document.getElementById('viz');
const ctx = canvas.getContext('2d');
const W = canvas.width, H = canvas.height;
const status = id('status');
const diag = id('diag');

function id(n){return document.getElementById(n)}

// Complex arithmetic
function C(re,im=0){return {re,im}}
function add(a,b){return C(a.re+b.re,a.im+b.im)}
function sub(a,b){return C(a.re-b.re,a.im-b.im)}
function mul(a,b){return C(a.re*b.re - a.im*b.im, a.re*b.im + a.im*b.re)}
function scale(a,s){return C(a.re*s,a.im*s)}
function abs(a){return Math.hypot(a.re,a.im)}
function conj(a){return C(a.re,-a.im)}

// Evaluate polynomial p(x) = x^4 + m x^3 + m^2 x^2 + m^3 x + m^4
function polyEval(x,m){
  // x is complex
  const x2 = mul(x,x), x3=mul(x2,x), x4=mul(x3,x);
  const m1 = m;
  const m2 = mul(m,m);
  const m3 = mul(m2,m);
  const m4 = mul(m2,m2);
  // x4 + m x3 + m2 x2 + m3 x + m4
  let s = add(x4, mul(m1,x3));
  s = add(s, mul(m2,x2));
  s = add(s, mul(m3,x));
  s = add(s, m4);
  return s;
}

// Durand-Kerner for quartic
function rootsDurandKerner(m, maxIter=200, tol=1e-12){
  // initial guesses on circle
  const n=4; const roots = [];
  const R = 1 + abs(m.re || 0) + abs(m.im || 0);
  for(let k=0;k<n;k++){
    const theta = 2*Math.PI*k/n;
    roots.push(C(R*Math.cos(theta), R*Math.sin(theta)));
  }
  for(let it=0; it<maxIter; it++){
    let moved=0;
    for(let i=0;i<n;i++){
      let xi = roots[i];
      let pxi = polyEval(xi,m);
      let denom = C(1,0);
      for(let j=0;j<n;j++) if(i!==j){
        denom = mul(denom, sub(xi, roots[j]));
      }
      // avoid division by tiny denom
      const denomAbs = abs(denom) + 1e-18;
      const delta = C(pxi.re/denom.re || pxi.re/denomAbs, pxi.im/denom.re || pxi.im/denomAbs);
      // Newton-like correction: xi = xi - p(xi)/prod(xi-xj)
      roots[i] = sub(xi, delta);
      moved = Math.max(moved, abs(delta));
    }
    if(moved < tol) break;
  }
  return roots;
}

// Line through two complex points returns {a,b,c} for ax+by+c=0 in real plane
function lineFromPoints(z1,z2){
  const x1=z1.re,y1=z1.im,x2=z2.re,y2=z2.im;
  const a = y1 - y2;
  const b = x2 - x1;
  const c = x1*y2 - x2*y1;
  return {a,b,c};
}
function intersectLines(L1,L2){
  const D = L1.a*L2.b - L2.a*L1.b;
  if(Math.abs(D) < 1e-12) return null;
  const x = (L1.b*L2.c - L2.b*L1.c)/D;
  const y = (L2.a*L1.c - L1.a*L2.c)/D;
  return C(x,y);
}
function collinearityDet(z1,z2,z3){
  const m = z1.re*(z2.im - z3.im) + z2.re*(z3.im - z1.im) + z3.re*(z1.im - z2.im);
  return m;
}

// Visualization helpers
function worldToCanvas(z){
  // center and scale
  const cx = W/2, cy = H/2;
  const scale = Math.min(W,H)/6;
  return {x:cx + z.re*scale, y:cy - z.im*scale};
}

function draw(){
  ctx.clearRect(0,0,W,H);
  // draw axes
  ctx.strokeStyle = 'rgba(255,255,255,0.06)'; ctx.lineWidth=1;
  ctx.beginPath(); ctx.moveTo(W/2,0); ctx.lineTo(W/2,H); ctx.moveTo(0,H/2); ctx.lineTo(W,H/2); ctx.stroke();

  // compute m
  const R = parseFloat(id('mag').value);
  const th = parseFloat(id('angle').value);
  const m = C(R*Math.cos(th), R*Math.sin(th));
  const roots = rootsDurandKerner(m);
  // map root indices to two base lines
  // L1: r0,r1,r2 ; L2: r1,r2,r3  (example overlapping choice)
  const r0=roots[0], r1=roots[1], r2=roots[2], r3=roots[3];

  // draw roots
  const points = [r0,r1,r2,r3];
  points.forEach((z,i)=>{
    const p = worldToCanvas(z);
    ctx.fillStyle = ['#ffd166','#ef476f','#06d6a0','#118ab2'][i%4];
    ctx.beginPath(); ctx.arc(p.x,p.y,6,0,Math.PI*2); ctx.fill();
    ctx.fillStyle='rgba(255,255,255,0.7)'; ctx.font='12px monospace';
    ctx.fillText('r'+i, p.x+8, p.y-6);
  });

  // chosen line points for Pappus-like construction
  const A1=r0, B1=r1, C1=r2;
  const A2=r1, B2=r2, C2=r3;

  // draw lines L1 and L2
  const L1 = lineFromPoints(A1,C1);
  const L2 = lineFromPoints(A2,C2);
  drawLine(L1,'rgba(110,231,183,0.25)');
  drawLine(L2,'rgba(92,122,255,0.14)');

  // Pappus pairing: intersections of (A1B2 with A2B1), (B1C2 with B2C1), (C1A2 with C2A1)
  const AB = lineFromPoints(A1,B2);
  const BA = lineFromPoints(A2,B1);
  const P = intersectLines(AB,BA);

  const BC = lineFromPoints(B1,C2);
  const CB = lineFromPoints(B2,C1);
  const Q = intersectLines(BC,CB);

  const CA = lineFromPoints(C1,A2);
  const AC = lineFromPoints(C2,A1);
  const Rpt = intersectLines(CA,AC);

  const triple = [P,Q,Rpt].filter(x=>x!==null);
  triple.forEach((z,i)=>{
    if(!z) return;
    const p = worldToCanvas(z);
    ctx.fillStyle = '#fffc9e';
    ctx.beginPath(); ctx.arc(p.x,p.y,5,0,Math.PI*2); ctx.fill();
    ctx.fillStyle='#071022'; ctx.font='11px monospace';
    ctx.fillText(['P','Q','R'][i], p.x+6, p.y+4);
  });

  // draw line through P,Q,R if collinear numerically
  if(triple.length===3){
    const det = collinearityDet(triple[0],triple[1],triple[2]);
    const eps = 1e-6;
    const ok = Math.abs(det) < eps;
    status.textContent = ok ? 'Collinear: YES' : 'Collinear: NO  (det='+det.toExponential(3)+')';
    if(ok){
      const Lp = lineFromPoints(triple[0], triple[1]);
      drawLine(Lp,'rgba(255,252,158,0.65)',2);
    } else {
      drawLine(lineFromPoints(triple[0], triple[1]), 'rgba(255,100,100,0.25)',1);
    }
    diag.textContent = [
      'm = '+m.re.toFixed(4)+' + '+m.im.toFixed(4)+'i',
      'roots:',
      ...roots.map((z,i)=>` r${i}: ${z.re.toFixed(5)} ${z.im>=0?'+':'-'} ${Math.abs(z.im).toFixed(5)}i`),
      'det(P,Q,R) = '+det.toExponential(6)
    ].join('\n');
  }
}

// draw infinite line from ax+by+c=0
function drawLine(L, color='rgba(255,255,255,0.08)', width=1){
  // compute two points on canvas
  const pts = [];
  if(Math.abs(L.b) > Math.abs(L.a)){
    // choose x extremes
    const x1 = -10, x2 = 10;
    const y1 = -(L.a*x1 + L.c)/L.b;
    const y2 = -(L.a*x2 + L.c)/L.b;
    pts.push(worldToCanvas(C(x1,y1)), worldToCanvas(C(x2,y2)));
  } else {
    const y1 = -10, y2 = 10;
    const x1 = -(L.b*y1 + L.c)/L.a;
    const x2 = -(L.b*y2 + L.c)/L.a;
    pts.push(worldToCanvas(C(x1,y1)), worldToCanvas(C(x2,y2)));
  }
  ctx.strokeStyle = color; ctx.lineWidth = width;
  ctx.beginPath(); ctx.moveTo(pts[0].x, pts[0].y); ctx.lineTo(pts[1].x, pts[1].y); ctx.stroke();
}

// animation / controls
let t0 = performance.now();
function stepLoop(){
  if(id('osc').checked){
    const w = parseFloat(id('omega').value);
    const now = performance.now();
    const dt = (now - t0)/1000;
    t0 = now;
    // advance angle by w * dt
    const angle = parseFloat(id('angle').value);
    id('angle').value = (angle + w*dt) % (2*Math.PI);
  }
  draw();
  requestAnimationFrame(stepLoop);
}

// wiring
id('step').addEventListener('click', ()=>{ draw(); });
id('reset').addEventListener('click', ()=>{ id('mag').value=1; id('angle').value=0; draw(); });

draw();
requestAnimationFrame(stepLoop);
