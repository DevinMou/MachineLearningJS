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
  let M = Array(rNum).fill(1);
  M = M.map((e, i) =>
    Array(cNum)
      .fill(null)
      .map((t, j) => {
        if (typeof fill === "function") {
          return fill(i, j);
        } else {
          return fill;
        }
      })
  );
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
    Q[0] = MPO(AT[0], R[0][0], MMath.div);
    // const min = ra < ca ? ra : ca;
    for (let j = 1; j < ca; j++) {
      let v = [...AT[j]];
      for (let i = 0; i < j; i++) {
        R[i][j] = MT(transpose(Q[i]), AT[j])[0][0];
        v = MPO(v, MPO(R[i][j], Q[i], MMath.time), MMath.sub);
      }
      R[j][j] = norm(v);
      Q[j] = MPO(v, R[j][j], MMath.div);
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
      Q[i] = MPO(V[i], R[i][i], MMath.div);
      for (let j = i + 1; j < ca; j++) {
        R[i][j] = MT(transpose(Q[i]), V[j])[0][0];
        V[j] = MPO(V[j], MPO(R[i][j], Q[i], MMath.time), MMath.sub);
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
    const na = norm(a);
    const b = createMatrix(rx, 1, (r, c) => (r === 0 && c === 0 ? na : 0));
    const t = MPO(a, b, MMath.sub);
    const th = MPO(
      createMatrix(rx, rx, (r, c) => (r === c ? 1 : 0)),
      MPO(MT(t, transpose(t)),2/norm(t,2,true), MMath.time),
      MMath.sub
    );
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

function Givens(A) {
  const [ra, ca] = getMRC(A);
  const g = (TX, a, b, i, T) => {
    const t1 = TX[i][a];
    const t2 = TX[i][b];
    let Q = null;
    if (t2 === 0) {
      return;
    } else if (t1 === 0) {
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

function ET (A) { // Elementary transformation
  const [ra, ca] = getMRC(A);
  if(ra!==null){
    const CA = copy(A)
    const min = ra < ca ? ra : ca;
    const tempA = []
    for(let i = 0 ; i<min;i++){
      const ni = CA.findIndex(e => !is0(e))
      if(ni!==-1){
        const tr = CA[ni].map(e=>e/CA[ni][i])
        CA.splice(ni,1)
        tempA.push(tr)
        CA.forEach(e=>{
          const coe = e[i]
          for(let j=i;j<ca;j++){
            e[j]-=coe*tr[j]
            is0(e[j])&&(e[j]=0);
          }
        })
      }
    }
    return [...tempA,...CA]
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

function vecI(v){ // vector identity
  const n = norm(v)
  return n ? v.map(e => e/n) : v
}

function getFeatVecByval (A, lam) {
  const HA = copy(A,(e,i,j)=>i===j?(e-lam):e)
  return SOE(ET(HA),true).slice(1).map(e=>vecI(e))
}
