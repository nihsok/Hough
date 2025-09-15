const omega = 2*Math.PI/86400 //earth rotation
const r_e = 6.37e6 //earth radius
const g = 9.81
const ev2ed = 4*r_e**2*omega**2/g
const n = 100 //truncation

function sigma2tau(){
  const sigma = document.getElementById('sigma').value
  document.getElementById('tau').value = 12/sigma
}

function tau2sigma(){
  const tau = document.getElementById('tau').value
  document.getElementById('sigma').value = 12/tau
}

function calc_eigen(){
  const s = parseFloat(document.getElementById('s').value)
  const sigma = parseFloat(document.getElementById('sigma').value)
  let M=[]
  let L=[]
  for (let r=s;r<=n;r++){ //r=s+i
    L[r-s] = Math.sqrt((r+s+1)*(r+s+2)*(r-s+1)*(r-s+2))/( (2*r+3)*Math.sqrt((2*r+1)*(2*r+5))*(s/sigma-(r+1)*(r+2)) )
    M[r-s] = sigma**2*(r*(r+1)-s/sigma)/( r**2*(r+1)**2 )
           + (r+2)**2*(r+s+1)*(r-s+1)/( (r+1)**2*(2*r+3)*(2*r+1)*(s/sigma-(r+1)*(r+2)) )
    if(s/sigma-r*(r-1) != 0) M[r-s]+= (r-1)**2*(r**2-s**2)/( r**2*(4*r**2-1)*(s/sigma-r*(r-1)) )
  }

  //symmetric mode
  const row = document.getElementsByName('F1_row')
  for (let i=0;i<4;i++){
    let j = i==0 ? 2+i : i //rowspan for the first row
    row[i].cells[j].innerText = M[i*2].toFixed(4)
    row[i].cells[j+1].innerText = L[i*2].toFixed(4)
    if(i<3) row[i+1].cells[i].innerText = L[i*2].toFixed(4)
  }

  const len_sym = s%2==0 ? (n-s+2)/2 : (n-s+1)/2
  const f1 = math.zeros(len_sym,len_sym)
  for (let i=0;i<len_sym;i++){
    f1.set([i,i],M[2*i])
    if (i+1<len_sym){
      f1.set([i,i+1],L[2*i])
      f1.set([i+1,i],L[2*i])
    }
  }

  const ans = math.eigs(f1)
  const lambda = ans.values._data
  const vector = ans.eigenvectors
  const x1 = document.getElementsByName('x1')
  for (i=1;i<10;i++){
    x1[0].cells[i+3].innerText =(lambda[len_sym-i]*ev2ed).toFixed(1)
    x1[1].cells[i+3].innerText = lambda[len_sym-i].toFixed(4)
    if (lambda[len_sym-i]<0){
      x1[0].cells[i+3].innerHTML = '<i>'+x1[0].cells[i+3].innerText+'</i>'
      x1[1].cells[i+3].innerHTML = '<i>'+x1[1].cells[i+3].innerText+'</i>'
    }
    x1[2].cells[i+3].innerText = vector[len_sym-i].vector._data[0].toFixed(4)
    x1[3].cells[i+1].innerText = vector[len_sym-i].vector._data[1].toFixed(4)
    x1[4].cells[i+1].innerText = vector[len_sym-i].vector._data[2].toFixed(4)
    x1[5].cells[i+1].innerText = vector[len_sym-i].vector._data[3].toFixed(4)
    x1[6].cells[i+1].innerText = '⋮'
  }
  //antisymmetric mode
}

//highlight column on hover
//plot