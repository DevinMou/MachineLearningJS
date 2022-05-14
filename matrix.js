/*
  A = LDU
  
*/
function GJ(A) {}
function LDU(A) {}

function norm(A, n = 2, nsqr) {
  let res = null;
  switch (n) {
    case 2:
      const [ra, ca] = getMRC(A);
      if (ra !== null && (ra === 1 || ca === 1)) {
        res = A.flat().reduce((a, b) => ((a += b ** 2), a), 0)
      }
      break;
  }
  return nsqr ? res : res**0.5;
}

function getMRC(A) {
  // get matrix number of row and column
  if (Array.isArray(A)) {
    const row = A.length;
    const isVector = !Array.isArray(A[0]);
    const col = isVector ? 1 : A[0].length;
    const err = A.find(
      (e) => !isVector ^ Array.isArray(e) || (!isVector && e.length !== col)
    );
    if (err) {
      return [null, null];
    } else {
      if (isVector) {
        A.forEach((e, i) => (A[i] = [e]));
      }
      return [row, col];
    }
  } else {
    return [null, null];
  }
}

function createMatrix(rNum, cNum, fill = null) {
  if (!rNum || !cNum) {
    return null;
  }
  const M = Array(rNum).fill(1).map(()=>Array(cNum).fill(1));
  const isFun = typeof fill === "function"
  M.forEach((r, i)=>{
    M[i] = r.map((e, j) => {
      if (isFun) {
        return fill(i, j, M);
      } else {
        return fill;
      }
    })
  })
  return M;
}

function transpose(A) {
  const [ra, ca] = getMRC(A);
  if (ra !== null) {
    let C = Array(ca).fill(1);
    C = C.map((e, i) => A.map((t) => t[i]));
    return C;
  } else {
    return null;
  }
}

function dot(A, B) {
  return A.reduce((a, b, c) => ((a += b * B[c]), a), 0);
}

function copy(A, fn) {
  return A.map((e,i)=>e.map((t,j)=>fn?fn(t,i,j):t))
}

function MT(A, B) {
  // matrix time
  const [ra, ca] = getMRC(A);
  const [rb, cb] = getMRC(B);
  if (ra !== null && rb !== null) {
    if (ca === rb) {
      let C = Array(ra).fill(1);
      const BT = transpose(B);
      C = C.map((e, i) => BT.map((b) => dot(b, A[i])));
      return C;
    }
  }
  return null;
}

function MPO(A, B, o) {
  const [ra, ca] = getMRC(A);
  const [rb, cb] = getMRC(B);
  if (typeof o === 'string') {
    o = MMath[MathSymbol[o]]
  }
  if (ra !== null && ra === rb && ca === cb) {
    const isVectorA = !Array.isArray(A[0]);
    const isVectorB = !Array.isArray(B[0]);
    let C = Array(ra).fill(1);
    C = C.map((a, i) =>
      Array(ca)
        .fill(1)
        .map((b, j) =>
          o(isVectorA ? A[i] : A[i][j], isVectorB ? B[i] : B[i][j])
        )
    );
    return C;
  } else if ((ra === null) ^ (rb === null)) {
    let C = Array(ra || rb).fill(1);
    const isMatrixA = ra !== null;
    const isVector = isMatrixA ? !Array.isArray(A[0]) : !Array.isArray(B[0]);
    C = C.map((a, i) =>
      Array(ca || cb)
        .fill(1)
        .map((b, j) =>
          o(
            isMatrixA ? (isVector ? A[i] : A[i][j]) : A,
            !isMatrixA ? (isVector ? B[i] : B[i][j]) : B
          )
        )
    );
    return C;
  } else {
    return null;
  }
}
const MathSymbol = {
  '*': 'time',
  '+': 'add',
  '/': 'div',
  '-': 'sub'
}
const MMath = {
  time: (a, b) => a * b,
  add: (a, b) => a + b,
  div: (a, b) => a / b,
  sub: (a, b) => a - b,
};

function GS(A) {
  const [ra, ca] = getMRC(A);
  if (ra !== null) {
    const AT = transpose(A);
    const R = createMatrix(ra, ca, 0);
    const Q = [];
    R[0][0] = norm(AT[0]);
    Q[0] = MPO(AT[0], R[0][0], '/');
    // const min = ra < ca ? ra : ca;
    for (let j = 1; j < ca; j++) {
      let v = [...AT[j]];
      for (let i = 0; i < j; i++) {
        R[i][j] = MT(transpose(Q[i]), AT[j])[0][0];
        v = MPO(v, MPO(R[i][j], Q[i], '*'), '-');
      }
      R[j][j] = norm(v);
      Q[j] = MPO(v, R[j][j], '/');
    }
    return [transpose(Q), R];
  } else {
    return [null, null];
  }
}

function MGS(A) {
  const [ra, ca] = getMRC(A);
  if (ra !== null) {
    const AT = transpose(A);
    const R = createMatrix(ra, ca, 0);
    const Q = [];
    const V = AT.map((e) => [...e]);
    const min = ra < ca ? ra : ca;
    for (let i = 0; i < min; i++) {
      R[i][i] = norm(V[i]);
      Q[i] = MPO(V[i], R[i][i], '/');
      for (let j = i + 1; j < ca; j++) {
        R[i][j] = MT(transpose(Q[i]), V[j])[0][0];
        V[j] = MPO(V[j], MPO(R[i][j], Q[i], '*'), '-');
      }
    }
    return [transpose(Q), R];
  } else {
    return [null, null];
  }
}

function splitM(A, [rs, re], [cs, ce]) {
  const [ra, ca] = getMRC(A);
  if (ra !== null) {
    re = re === undefined ? ra : re < 0 ? ra + re : re;
    ce = ce === undefined ? ca : ce < 0 ? ca + ce : ce;
    return A.slice(rs, re).map((e) => e.slice(cs, ce));
  } else {
    return null;
  }
}

function concatM(...As) { // [A B] (A,B,false)
  const last = As[As.length-1]
  const typeLast = typeof last
  const mode = typeLast === "boolean" ?  last ? 2 : 1 : 0
  mode && As.splice(As.length-1,1)
  const rcs = As.map(a=>getMRC(a))
  const equalR = !rcs.find(e => e[0]!==rcs[0][0])
  const equalC = !rcs.find(e => e[1]!==rcs[0][1])
  function sum (A) {
    return A.reduce((a,b)=>(a+=b,a))
  }
  if(mode ? mode===1 && equalR || mode===2 && equalC : equalR||equalC){
    const rcsT = transpose(rcs)
    const d = mode ? mode===1 : equalR
    return createMatrix(d? rcs[0][0] : sum(rcsT[0]),d? sum(rcsT[1]):rcs[0][1],(r,c)=>{
      if(d){
        const i = rcsT[1].findIndex(e=>c<e?true:(c-=e,false))
        return As[i][r][c]
      }else{
        const i = rcsT[0].findIndex(e=>r<e?true:(r-=e,false))
        return As[i][r][c]
      }
    })
  }else{
    return null
  }
}

function Hh(v) {
  const len = v.length
  if (len===1){
    return [v[0],[[1]]]
  }
  const nv = norm(v)
  if (len === 2) {
    const s = v[0] / nv
    const c = v[1] / nv
    return [[nv,0],[[s,c],[c,-s]]]
  }
  const w = createMatrix(len,1,(r,c)=>r===0&&c===0?nv:0)
  const t = MPO(v,w, '-')
  const P = MPO(createMatrix(len,len,(r,c)=>r===c?1:0),MPO(MT(t,transpose(t)),2/norm(t,2,true),'*'),'-')
  return [w, P]
}

function HR(A, getQ) { // Hessenburg Reduction
  const [ra, ca] = getMRC(A)
  let TA = copy(A)
  let P = null
  for(let i=0;i<ca-1;i++){
    const v = TA.slice(1+i-ca).map(e=>e[i])
    const [w,p] = Hh(v)
    const tP = createMatrix(ra,ca,(r,c)=>r>i&&c>i?(p[r-1-i][c-1-i]):r===c?1:0)
    if (getQ) {
      P = P ? MT(P,tP) : tP
    } 
    const tp = transpose(p)
    const X0 = splitM(TA,[0,i+1],[0,i+1])
    const X1 = MT(p,splitM(TA,[i+1],[0,i+1]))
    const X2 = MT(MT(p,splitM(TA,[i+1],[i+1])),tp)
    const X3 = MT(splitM(TA,[0,i+1],[i+1]),tp)
    TA = concatM(concatM(X0,X3),concatM(X1,X2),true)
  }
  return [P,TA]
}

function HH(A) {
  // Householder
  const [ra, ca] = getMRC(A);
  const min = ra < ca ? ra : ca;
  let H = null;
  let AA = A;
  for (let level = 0; level < min; level++) {
    const X = splitM(AA, [level], [level]);
    const rx = X.length;
    const XT = transpose(X);
    const a = XT[0];
    const [b,th] = Hh(a)
    const [rh, ch] = getMRC(th);
    const h = createMatrix(ra, ra, (r, c) => {
      if (r >= ra - rh && c >= ra - rh) {
        return th[r - ra + rh][c - ra + rh];
      } else if (r === c) {
        return 1;
      } else {
        return 0;
      }
    });
    if (H) {
      H = MT(H, h);
    } else {
      H = h;
    }
    AA = MT(h, AA);
  }
  return [H, AA];
}

function givens(x, y) {
  const z = Math.sqrt(Math.pow(x, 2) + Math.pow(y, 2))
  return z ? [x / z, y / z] : []
}

function Givens(A) {
  const [ra, ca] = getMRC(A);
  const g = (TX, a, b, i, T) => {
    const t1 = TX[i][a];
    const t2 = TX[i][b];
    let Q = null;
    if (is0(t2)) {
      return;
    } else if (is0(t1)) {
      TX.forEach((e) => {
        [e[a], e[b]] = [e[b], e[a]];
      });
      Q = createMatrix(ra, ra, (r, c) => {
        if ((r === a && c === a) || (r === b && c === b)) {
          return 0;
        } else if (r === a && c === b) {
          return 1;
        } else if (r === b && c === a) {
          return -1;
        } else if (r === c) {
          return 1;
        } else {
          return 0;
        }
      });
    } else {
      const n = Math.sqrt(Math.pow(t1, 2) + Math.pow(t2, 2));
      const cos = t1 / n;
      const sin = t2 / n;
      TX.forEach((e) => {
        const ea = e[a];
        const eb = e[b];
        e[a] = cos * ea + sin * eb;
        e[b] = -sin * ea + cos * eb;
      });
      Q = createMatrix(ra, ra, (r, c) => {
        if ((r === a && c === a) || (r === b && c === b)) {
          return cos;
        } else if (r === a && c === b) {
          return sin;
        } else if (r === b && c === a) {
          return -sin;
        } else if (r === c) {
          return 1;
        } else {
          return 0;
        }
      });
    }
    if (!T.Q) {
      T.Q = Q;
    } else {
      T.Q = MT(Q, T.Q);
    }
  };
  if (ra !== null) {
    const AT = transpose(A);
    let T = { Q: null };
    for (let i = 0; i < ca; i++) {
      for (let j = ra - 2; j >= i; j--) {
        g(AT, j, j + 1, i, T);
      }
    }
    return [transpose(T.Q), transpose(AT)];
  }
}

function is0 (n) { // equalZero
  return n < 0 ? n > -1e-12 : n < 1e-12
}

function ET (A, reduce, augmented) { // Elementary transformation
  const [ra, ca] = getMRC(A);
  if(ra!==null){
    const [ram,cam] = augmented ? getMRC(augmented) : [0,0]
    const CA = ram ? concatM(A,augmented) : copy(A)
    const min = ra < ca ? ra : ca;
    const tempA = []
    for(let i = 0 ; i<min;i++){
      const ni = CA.findIndex((e, j) => j<ca && !is0(e[i]))
      if(ni!==-1){
        const tr = CA[ni].map(e=>e/CA[ni][i])
        CA.splice(ni,1)
        tempA.push(tr)
        CA.forEach(e=>{
          const coe = e[i]
          for(let j=i;j<ca+cam;j++){
            e[j]-=coe*tr[j]
            is0(e[j])&&(e[j]=0);
          }
        })
      }
    }
    if(reduce){
      for(let i=1;i<tempA.length;i++){
        const item = tempA[i]
        const ni = item.findIndex(e => !is0(e))
        for(let j=0;j<i;j++){
          const n = tempA[j][ni]/item[ni]
          tempA[j] = tempA[j].map((e,k)=>e-n*item[k])
        }
      }
    }
    if(CA.find(e=>e.find((t,i)=>i<ca && !is0(t)))){
      return null
    }
    return augmented ? [[splitM(tempA,[0,ra],[0,ca]),splitM(CA,[0,ra],[0,ca])],[splitM(tempA,[0,ra],[ca]),splitM(CA,[0,ra],[ca])]] : [tempA,CA]
  }
  return null
}

function SOE (E, h) { // solution of equation, h: Homogeneous linear
  E = copy(E)
  if(h){
    E.forEach(e=>e.push(0))
  }
  const len = E[0].length
  const tA = Array(len-1).fill(null)
  E.forEach(e => {
    const i = e.findIndex(t=>t===1)
    if(i>-1){
      tA[len-2-i] = [e[len-1],e.slice(i+1,-1).reverse().map(t=>-t)]
    }
  })
  const RA = [[]]
  tA.forEach((e, i)=>{
    if(e){
      if(e.length===1&&!RA[0][0].length){
        RA[0].push(e[0])
      }else{
        MT(RA,transpose([e[1]])).forEach((a,b)=>RA[b].push(b?a[0]:(a[0]+e[0])))
      }
    }else{
      RA.forEach(a=>a.push(0))
      RA.push([...Array(i).fill(0),1])
    }
  })
  return RA.map(e => e.reverse())
}

function AXB (E, hasB) { // resolve Ax=b
  E = copy(E)
  const len = E[0].length + (hasB ? 1 : 0)
  let r = len - 1
  const RI = Object.keys(Array(len-1).fill()).map(e=>+e)
  E.forEach(e=>{
    const i = e.findIndex(t => !is0(t))
    if(i>-1){
      RI.splice(RI.indexOf(i),1)
      r-=1
    }
    hasB&&e.push(0)
  })
  RI.forEach(e => {
    const t = Array(len).fill(0)
    t[e] = 1
    E.splice(e,0,t)
  })
  const RA = transpose(E.slice(0,len-1)).filter((e, i) => RI.includes(i)||i===len-1)
  return [RA[RA.length-1],...RA.slice(0,-1).map(e=>e.map(t=>-t))]
}

function vecI(v){ // vector identity
  const n = norm(v)
  return n ? v.map(e => e/n) : v
}

function getFeatVecByval (A, lam) {
  const HA = copy(A,(e,i,j)=>i===j?(e-lam):e)
  return SOE(ET(HA),true).slice(1).map(e=>vecI(e))
}

function fixed(A,n=4,tn) {
  return A.map(e=>e.map(t=>(s=t.toFixed(n),tn?+s:s)))
}

function ws (A, n=1) {
  let H = HR(A)[1]
  const ra = getMRC(H)[0] - 1
  while(n--){
    const s = H[ra][ra] + H[ra-1][ra-1]
    const t = H[ra][ra]*H[ra-1][ra-1] - H[ra][ra-1]*H[ra-1][ra]
    const M = MPO(MPO(MT(H,H),MPO(s,H,'*'),'-'),createMatrix(ra+1,ra+1,(r,c)=>r===c?t:0),'+')
    const [Q,R] = Givens(M)
    H = MT(MT(transpose(Q),H),Q)
  }
  return H
}

function fs (A, n=1) {
  let H = HR(A)[1]
  const ra = getMRC(H)[0] - 1
  while(n--){
    const s = H[ra][ra] + H[ra-1][ra-1]
    const t = H[ra][ra]*H[ra-1][ra-1] - H[ra][ra-1]*H[ra-1][ra]
    const x = H[0][0]**2 + H[0][1]*H[1][0] - s*H[0][0] + t
    const y = H[1][0]*(H[0][0]+H[1][1]-s)
    const z = H[1][0]*H[2][1]
    const p = Hh([x,y,z].concat(Array(ra-2).fill(0)))[1]
    const P = HR(MT(MT(transpose(p),H),p), true)[0]
    const Q = MT(p,P)
    H = MT(MT(transpose(Q),H),Q)
  }
  return H
}

function replaceM (A, B, [rs,re],[cs,ce]) {
  for(let r=rs;r<re;r++){
    for(let c=cs;c<ce;c++){
      A[r][c] = B[r-rs][c-cs]
    }
  }
}

function sfs (A, tol=1) {
  let Hb = HR(A)[1]
  const H=(a,b)=>Hb[a-1][b-1]
  const n = getMRC(Hb)[0]
  let p = n
  const e = 1e-15
  const log = []
  let time = 0
  while(tol--&&p>2){
    let q = p -1
    let s = H(q,q) + H(p,p)
    let t = H(q,q)*H(p,p) - H(p,q)*H(q,p)
    let x = H(1,1)**2 + H(1,2)*H(2,1) - s*H(1,1) + t
    let y = H(2,1)*(H(1,1)+H(2,2)-s)
    let z = H(2,1)*H(3,2)
    for(let k=0;k<p-2;k++){
      let P = Hh([x,y,z])[1]
      let tP = transpose(P)
      let r = Math.max(1,k)
      replaceM(Hb,MT(tP,splitM(Hb,[k,k+3],[r-1,n])),[k,k+3],[r-1,n])
      r = Math.min(k+4,p)
      replaceM(Hb,MT(splitM(Hb,[0,r],[k,k+3]),P),[0,r],[k,k+3])
      x = H(k+2,k+1)
      y = H(k+3,k+1)
      if (k < p-3) {
        z = H(k+4,k+1)
      }
    }
    let P = Hh([x,y])[1]
    let tP = transpose(P)
    replaceM(Hb,MT(tP,splitM(Hb,[q-1,p],[p-3,n])),[q-1,p],[p-3,n])
    replaceM(Hb,MT(splitM(Hb,[0,p],[p-2,p]),P),[0,p],[p-2,p])
    time++
    log.push([time, p, H(p,p-1),H(p-1,p-2)])
    if(Math.abs(H(p,q))<e*(Math.abs(H(q,q))+Math.abs(H(p,p)))){
      Hb[p-1][q-1] = 0
      p = p - 1
      q = p - 1
    } else if (Math.abs(H(p-1,q-1))<e*(Math.abs(H(q-1,q-1))+Math.abs(H(q,q)))) {
      Hb[p-2][q-2] = 0
      p = p - 2
      q = p - 1
    }
  }
  console.log(tol)
  return [Hb, Hb.map((e,i)=>e[i])]
}

function stqr (oA, tol = 1) {
  const A = HR(oA)[1]
  const n = getMRC(A)[0]
  let m = n
  const e = 1e-15
  const Q = createMatrix(n, n, (r,c) => r===c ? 1 : 0)
  const a = (i, v) => v === undefined ? A[i-1][i-1] : (A[i-1][i-1] = v, true)
  const b = (i, v) => v === undefined ? A[i-1][i-2] : (A[i-1][i-2] = v, A[i-2][i-1] = v, true)
  let d,s,x,y,gc,gs,w,z
  while(tol-- && m > 1) {
    d = (a(m-1) - a(m)) / 2
    if (d === 0) {
      s = a(m) - Math.abs(b(m))
    } else {
      const b2 = Math.pow(b(m),2)
      s = a(m) - b2 / (d + Math.sign(d)*Math.sqrt(b2+Math.pow(d,2)))
    }
    x = a(1) - s
    y = b(2)
    for(let k = 1; k <= m - 1; k++){
      if (m > 2) {
        [gc, gs] = givens(x, y)
      } else {
        const t = Math.pow((a(1) - a(2)) / b(2), 2)
        const c2 = (4 + t - Math.sqrt(t*(t+4))) / 2
        gc = Math.sqrt(c2)
        gs = Math.sqrt(1 - c2)
      }
      w = gc*x - gs*y
      d = a(k) - a(k+1)
      z = (2*gc*b(k+1) + d*gs)*gs
      a(k, a(k) - z)
      a(k + 1, a(k + 1) + z)
      b(k + 1, d*gc*gs + (Math.pow(gc, 2) - Math.pow(gs, 2))*b(k + 1))
      x = b(k + 1)
      if (k > 1) {
        b(k, w)
      }
      if (k < m - 1) {
        y = -gs*b(k + 2)
        b(k + 2, gc*b(k + 2))
      }
      replaceM(Q, MT(splitM(Q,[0,n],[k-1,k+1]),[[gc, gs],[-gs,gc]]),[0,n],[k-1,k+1])
    }
    const bm = Math.abs(b(m))
    const am = Math.abs(a(m-1)) + Math.abs(a(m))
    if (bm < e*am) {
      m = m - 1
    }
  }
  console.log(tol)
  return [A, Q]
}
function stqrc (oA, tol = 1) {
  const A = HR(oA)[1]
  const n = getMRC(A)[0]
  const Q = createMatrix(n, n, (r,c) => r===c ? 1 : 0)
  const e = 1e-15
  let d,s,x,y
  for(let m=n;m>=2;m--){
    while(1) {
      const cs = Array(n-1).fill(0)
      const ss = Array(n-1).fill(0)
      d = (A[m-2][m-2] - A[m-1][m-1]) / 2
      const b2 = Math.pow(A[m-1][m-2],2)
      s = A[m-1][m-1] - Math.sign(d) * b2 / (Math.abs(d) + Math.sqrt(b2+Math.pow(d,2)))
      for(let i=0; i < m; i++){
        A[i][i] -=s
      }
      for(let i=0; i < m-1; i++){
        x = A[i][i]
        y = A[i+1][i]
        const [_c, _s] = givens(x, y)
        cs[i] = _c
        ss[i] = _s
        for(let j=i;j<i+2 && j<m;j++){
          x = A[i][j]
          y = A[i+1][j];
          A[i][j] = cs[i]*x+ss[i]*y
          A[i+1][j] = -ss[i]*x+cs[i]*y
        }
      }
      for(let j=1;j<m;j++){
        for(let i=Math.max(0,j-2);i<=j;i++){
          x=A[i][j-1]
          y=A[i][j]
          A[i][j-1]=cs[j-1]*x+ss[j-1]*y
          A[i][j]=-ss[j-1]*x+cs[j-1]*y
        }
      }
      for(let i=0;i<m;i++){
        A[i][i]+=s
      }
      for(let j=1;j<m;j++){
        for(let i=0;i<n;i++){
          x=Q[i][j-1]
          y=Q[i][j]
          Q[i][j-1]=cs[j-1]*x+ss[j-1]*y
          Q[i][j]=-ss[j-1]*x+cs[j-1]*y
        }
      }
      const bm = Math.abs(A[m-1][m-2])
      const am = Math.abs(A[m-2][m-2]) + Math.abs(A[m-1][m-1])
      if((bm < e*am)){
        break
      }
    }
  }
  console.log(tol)
  return [A, Q]
}

function deflation(A) {
  const tempT = []
  const tempI = []
  const tempR = []
  const len = A.length
  for(let i = 0; i < len-1;i++){
    if (is0(A[i+1][i])) {
      tempT.push(A[i][i])
      tempI.push(i)
    } else {
      tempR.push(i)
    }
  }
  return [tempR.map(i => tempR.map(j => A[i][j])), tempT, tempI]
}

function strMM(A, B, toStr) {
  const splStr = (str) => {
    const m = str.match(/(?:\+|-)?\w+(?:\*(?:\+|-)?\w+)*/g)
    if (m && m.length) {
      return m.map(e => {
        let s = true
        let n = 1
        const res = e.split('*').map(t => {
          const m$ = /^(\+|-)?(\w+)/.exec(t)
          s = m$[1] === '-' ? !s : s
          if (/^\d+$/.test(m$[2])){
            n*=m$[2]
            return null
          } else {
            return m$[2]
          }
        }).filter(t => t!==null)
        return [s ? n : -n, res]
      })
    } else {
      return null
    }
  }
  const mul = (a,r,b,c) => {
    const lr = a[r]
    const rr = b.map((e)=>e[c])
    const res = lr.map((a,i)=>{
      const b = rr[i]
      return a.map(([na,resa])=>{
        return b.map(([nb,resb])=>{
          const arr = [...resa,...resb].sort()
          if (arr.find(e=>/^0+$/.test(e)) || na*nb===0){
            return 0
          } else {
            return [na*nb,arr]
          }
        })
      })
    }).flat(2).filter(e => e)
    const dist = {}
    const arr = []
    res.forEach((e,i)=>{
      const k = e[1].join('*')
      dist[k] = dist[k] || []
      !dist[k].length && arr.push(k)
      dist[k].push(i)
    })
    return arr.map(k => {
      return [dist[k].reduce((a,b)=>(a+=res[b][0],a),0),res[dist[k][0]][1]]
    })
  }
  const PA = A.map(e => e.map(t => splStr(t)))
  const PB = B.map(e => e.map(t => splStr(t)))
  const C = createMatrix(getMRC(A)[0],getMRC(B)[1],(r,c)=>{
    const res = mul(PA,r,PB,c)
    if (toStr) {
      return res.map((e,i)=>{
        const abs = Math.abs(e[0])
        return `${e[0]>0?'+':'-'}${abs===1?'':(abs+(e[1].length?'*':''))}${e[1].join('*')}`
      }).join('').replace(/^\+/,'') || '0'
    } else {
      return res
    }
  })
  return C
} 