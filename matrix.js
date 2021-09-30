/*
  A = LDU
  
*/
function GJ(A) {
  
}
function LDU(A) {
  

}

function norm(A, n=2){
  let res = null
  switch(n){
    case 2:
      const [ra,ca] = getMRC(A)
      if(ra!==null&&(ra===1||ca===1)){
        res = A.flat().reduce((a,b)=>(a+=b**2,a),0)**0.5
      }
      break
  }
  return res
}

function getMRC(A) { // get matrix number of row and column
  if(Array.isArray(A)){
    const row = A.length
    const isVector = !Array.isArray(A[0])
    const col = isVector? 1 : A[0].length
    const err = A.find(e => !isVector ^ Array.isArray(e) || !isVector && e.length !== col)
    if(err){
      return [null,null]
    }else{
      if(isVector){
        A.forEach((e,i)=>A[i] = [e])
      }
      return [row, col]
    }
  }else{
    return [null,null]
  }
}

function createMatrix(rNum,cNum,fill=null){
  if(!rNum||!cNum){
    return null
  }
  let M = Array(rNum).fill(1)
  M = M.map((e, i) => Array(cNum).fill(null).map((t, j)=>{
    if(typeof fill === 'function'){
      return fill(i,j)
    }else{
      return fill
    }
  }))
  return M
}

function transpose(A) {
  const [ra,ca] = getMRC(A)
  if(ra!==null){
    let C = Array(ca).fill(1)
    C = C.map((e, i)=>A.map(t=>t[i]))
    return C
  }else{
    return null
  }
}

function dot(A, B) {
  return A.reduce((a,b,c) => (a+=b*B[c],a),0)
}

function MT(A, B) { // matrix time
  const [ra,ca] = getMRC(A)
  const [rb,cb] = getMRC(B)
  if(ra!==null&&rb!==null){
    if(ca===rb){
      let C = Array(ra).fill(1)
      const BT = transpose(B)
      C = C.map((e, i)=>BT.map(b=>dot(b, A[i])))
      return C
    }
  }
  return null
}

function MPO(A,B,o) {
  const [ra,ca] = getMRC(A)
  const [rb,cb] = getMRC(B)
  if(ra!==null&&ra===rb&&ca===cb){
    const isVectorA = !Array.isArray(A[0])
    const isVectorB = !Array.isArray(B[0])
    let C = Array(ra).fill(1)
    C = C.map((a,i)=>Array(ca).fill(1).map((b,j)=>o(isVectorA ? A[i] : A[i][j],isVectorB ? B[i]:B[i][j])))
    return C
  }else if(ra===null^rb===null){
    let C = Array(ra||rb).fill(1)
    const isMatrixA = ra!==null
    const isVector = isMatrixA ? !Array.isArray(A[0]) : !Array.isArray(B[0])
    C = C.map((a,i)=>Array(ca||cb).fill(1).map((b,j)=>o(isMatrixA ? isVector ? A[i] : A[i][j] : A,!isMatrixA ? isVector ? B[i] : B[i][j] : B)))
    return C
  }else{
    return null
  }
}

const MMath = {
  time: (a,b)=>a*b,
  add: (a,b)=>a+b,
  div: (a,b)=>a/b,
  sub: (a,b)=>a-b
}

function GS(A) {
  const [ra,ca] = getMRC(A)
  if(ra!==null){
    const AT = transpose(A)
    const R = createMatrix(ca,ca,0)
    const Q = []
    R[0][0] = norm(AT[0])
    Q[0] = MPO(AT[0],R[0][0],MMath.div)
    for(let j=1;j<ca;j++){
      let v = [...AT[j]]
      for(let i=0;i<j;i++){
        R[i][j] = MT(transpose(Q[i]),AT[j])[0][0]
        v = MPO(v,MPO(R[i][j],Q[i],MMath.time),MMath.sub)
      }
      R[j][j] = norm(v)
      Q[j] = MPO(v,R[j][j],MMath.div)
    }
    return [transpose(Q),R]
  }else{
    return [null,null]
  }
}

function MGS(A) {
  const [ra,ca] = getMRC(A)
  if(ra!==null){
    const AT = transpose(A)
    const R = createMatrix(ca,ca,0)
    const Q = []
    const V = AT.map(e=>[...e])
    for(let i = 0;i<ca;i++){
      R[i][i] = norm(V[i])
      Q[i] = MPO(V[i],R[i][i],MMath.div)
      for(let j = i+1;j<ca;j++){
        R[i][j] = MT(transpose(Q[i]),V[j])[0][0]
        V[j] = MPO(V[j],MPO(R[i][j],Q[i],MMath.time),MMath.sub)
      }
    }
    return [transpose(Q),R]
  }else{
    return [null,null]
  }
}

function splitM (A, [rs,re], [cs,ce]){
  const [ra, ca] = getMRC(A)
  if(ra!==null){
    re = re === undefined ? ca : re < 0 ? ca + re : re
    ce = ce === undefined ? ra : ce < 0 ? ra + ce : ce
    return A.slice(rs,re).map(e=>e.slice(cs,ce))
  }else{
    return null
  }
}

function HH(A) { // Householder
  const [ra, ca] = getMRC(A)
  const min = ra < ca ? ra : ca
  let H = null
  function hh (AA, level=0) {
    if(level>min-2){
      return AA
    }
    const X = splitM(AA,[level],[level])
    const rx = X.length
    const XT = transpose(X)
    const a = XT[0]
    const na = norm(a)
    const b = createMatrix(rx,1,(r,c)=>r===0&&c===0?na:0)
    const t = MPO(a,b,MMath.sub)
    const w = MPO(t,norm(t),MMath.div)
    const th = MPO(createMatrix(rx,rx,(r,c)=>r===c?1:0), MPO(2, MT(w,transpose(w)), MMath.time),MMath.sub)
    const [rh, ch] = getMRC(th)
    const h = createMatrix(ra, ca,(r,c)=>{
      if((r>=ra - rh)&&(c>=ca - ch)){
        return th[r-ra+rh][c-ca+ch]
      }else if(r===c){
        return 1
      }else{
        return 0
      }
    })
    if(H){
      H = MT(H,h)
    }else{
      H = h
    }
    return hh(MT(h,AA),level+1)
  }
  const R = hh(A)
  return [H,R]
}