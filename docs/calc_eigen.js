const omega = 2*Math.PI/86400 //earth rotation
const r_e = 6.37e6 //earth radius
const g = 9.81
const ev2ed = 4*r_e**2*omega**2/g
const rmax = 100 //truncation

function sigma2tau(){
  const sigma = document.getElementById('sigma').value
  document.getElementById('tau').value = 12/sigma
}
function tau2sigma(){
  const tau = document.getElementById('tau').value
  document.getElementById('sigma').value = 12/tau
}

function calc_eigen(){
  const s = parseInt(document.getElementById('s').value)
  const sigma = parseFloat(document.getElementById('sigma').value)
  let M=[]
  let L=[]
  for (let r=s;r<=rmax;r++){ //r=s+i
    L[r-s] = Math.sqrt((r+s+1)*(r+s+2)*(r-s+1)*(r-s+2))
           /( (2*r+3)*Math.sqrt((2*r+1)*(2*r+5))*(s/sigma-(r+1)*(r+2)) )
    M[r-s] = sigma**2*(r*(r+1)-s/sigma)/( r**2*(r+1)**2 )
           + (r+2)**2*(r+s+1)*(r-s+1)/( (r+1)**2*(2*r+3)*(2*r+1)*(s/sigma-(r+1)*(r+2)) )
    if(s/sigma-r*(r-1) != 0) M[r-s]+= (r-1)**2*(r**2-s**2)/( r**2*(4*r**2-1)*(s/sigma-r*(r-1)) )
  }

  //symmetric mode
  const row_F1 = document.getElementsByName('F1_row')
  const len_sym = s%2==0 ? (rmax-s+2)/2 : (rmax-s+1)/2 //this should be constant, for the precision
  const f1 = math.zeros(len_sym,len_sym)
  for (let i=0;i<len_sym;i++){
    f1.set([i,i],M[2*i])
    if (i+1<len_sym){
      f1.set([i,i+1],L[2*i])
      f1.set([i+1,i],L[2*i])
    }
    if (i<5){
      const j = i===0 ? i+2 : i //rowspan for the first row
      row_F1[i].cells[j].innerText = M[i*2].toFixed(3)
      if(i===4) continue
      row_F1[i].cells[j+1].innerText = L[i*2].toFixed(3)
      row_F1[i+1].cells[i].innerText = L[i*2].toFixed(3)
    }
  }

  let ans = math.eigs(f1)
  const lambda_x1 = ans.values._data
  const vector_x1 = ans.eigenvectors
  const table_x1 = document.getElementById('X1')
  const table_x2 = document.getElementById('X2')
  let a_sym = []
  let n_p = s-2
  let n_n = s
  for (let i=1;i<=12;i++){
    table_x1.rows[0].cells[i+2].innerHTML = 'Θ<sub>'+(lambda_x1[len_sym-i]>0 ? n_p+=2 : n_n-=2)+'</sub><sup>'+s+','+sigma+'</sup>'
    table_x1.rows[1].cells[i+2].innerText = (4*Math.PI*7/Math.sqrt(8000/(lambda_x1[len_sym-i]*ev2ed)-1)).toFixed(1)
    table_x1.rows[2].cells[i+2].innerText =(lambda_x1[len_sym-i]*ev2ed).toFixed(1)
    table_x1.rows[3].cells[i+2].innerText = lambda_x1[len_sym-i].toFixed(4)
    if (lambda_x1[len_sym-i]<0){
      table_x1.rows[2].cells[i+2].innerHTML = '<i>'+table_x1.rows[2].cells[i+2].innerText+'</i>'
      table_x1.rows[3].cells[i+2].innerHTML = '<i>'+table_x1.rows[3].cells[i+2].innerText+'</i>'
    }
    table_x1.rows[4].cells[i+2].innerText = vector_x1[len_sym-i].vector._data[0].toExponential(2)
    table_x1.rows[5].cells[i  ].innerText = vector_x1[len_sym-i].vector._data[1].toExponential(2)
    table_x1.rows[6].cells[i  ].innerText = vector_x1[len_sym-i].vector._data[2].toExponential(2)
    table_x1.rows[7].cells[i  ].innerText = vector_x1[len_sym-i].vector._data[3].toExponential(2)
    table_x1.rows[8].cells[i  ].innerText = vector_x1[len_sym-i].vector._data[4].toExponential(2)
    table_x1.rows[9].cells[i  ].innerText = '⋮'

    a_sym[i-1] = vector_x1[len_sym-i].vector._data
  }

  let currentindex = null
  table_x1.addEventListener("mouseover", e => {
    const cell = e.target.closest("td")
    if(!cell) return
    const table = e.target.closest('table')
    const colIndex = cell.cellIndex
    const rowIndex = cell.parentNode.rowIndex
    if(colIndex === 15) return
    if(rowIndex<5 && colIndex<3 || rowIndex>4 && colIndex<1) return
    if(colIndex-(rowIndex>4 ? 1:3) === currentindex) return
    table.querySelectorAll("td").forEach(cell=>{
      cell.classList.remove("highlight")
    })
    table_x2.querySelectorAll("td").forEach(cell=>{
      cell.classList.remove("highlight")
    })
    Array.from(table.rows).forEach((row,i)=>{
      row.cells[colIndex-(i>4 ? 2:0)+(rowIndex>4 ? 2:0)].classList.add("highlight")
    })
    currentindex = colIndex-(rowIndex>4 ? 1:3)
    plot_hough(s,sigma,0,a_sym[currentindex])
  })
  // table_x1.addEventListener('mouseleave', e=>{
  //   currentindex = null
  //   e.target.closest('table').querySelectorAll("td").forEach(cell=>{
  //     cell.classList.remove("highlight")
  //   })
  // })


  //antisymmetric mode //the other is not renewed problem
  const row_F2 = document.getElementsByName('F2_row')
  const len_asym = s%2==0 ? (rmax-s+2)/2-1 : (rmax-s+1)/2-1//?
  const f2 = math.zeros(len_asym,len_asym)
  for (let i=0;i<len_asym;i++){
    f2.set([i,i],M[2*i+1])
    if (i+1<len_asym){
      f2.set([i,i+1],L[2*i+1])
      f2.set([i+1,i],L[2*i+1])
    }
    if (i<5){
      const j = i===0 ? i+2 : i //rowspan for the first row
      row_F2[i].cells[j].innerText = M[i*2+1].toFixed(3)
      if(i===4) continue
      row_F2[i].cells[j+1].innerText = L[i*2+1].toFixed(3)
      row_F2[i+1].cells[i].innerText = L[i*2+1].toFixed(3)
    }
  }

  ans = math.eigs(f2)
  const lambda_x2 = ans.values._data
  const vector_x2 = ans.eigenvectors

  let a_asym = []
  let n_p_a = s-1
  let n_n_a = s-1
  for (let i=1;i<=12;i++){
    table_x2.rows[0].cells[i+2].innerHTML = 'Θ<sub>'+(lambda_x2[len_asym-i]>0 ? n_p_a+=2 : n_n_a-=2)+'</sub><sup>'+s+','+sigma+'</sup>'
    table_x2.rows[1].cells[i+2].innerText = (4*Math.PI*7/Math.sqrt(8000/(lambda_x2[len_asym-i]*ev2ed)-1)).toFixed(1)
    table_x2.rows[2].cells[i+2].innerText =(lambda_x2[len_asym-i]*ev2ed).toFixed(1)
    table_x2.rows[3].cells[i+2].innerText = lambda_x2[len_asym-i].toFixed(4)
    if (lambda_x2[len_asym-i]<0){
      table_x2.rows[2].cells[i+2].innerHTML = '<i>'+table_x2.rows[2].cells[i+2].innerText+'</i>'
      table_x2.rows[3].cells[i+2].innerHTML = '<i>'+table_x2.rows[3].cells[i+2].innerText+'</i>'
    }
    table_x2.rows[4].cells[i+2].innerText = vector_x2[len_asym-i].vector._data[0].toExponential(2)
    table_x2.rows[5].cells[i  ].innerText = vector_x2[len_asym-i].vector._data[1].toExponential(2)
    table_x2.rows[6].cells[i  ].innerText = vector_x2[len_asym-i].vector._data[2].toExponential(2)
    table_x2.rows[7].cells[i  ].innerText = vector_x2[len_asym-i].vector._data[3].toExponential(2)
    table_x2.rows[8].cells[i  ].innerText = vector_x2[len_asym-i].vector._data[4].toExponential(2)
    table_x2.rows[9].cells[i  ].innerText = '⋮'

    a_asym[i-1] = vector_x2[len_asym-i].vector._data
  }

  currentindex = null
  table_x2.addEventListener("mouseover", e => {
    const cell = e.target.closest("td")
    if(!cell) return
    const table = e.target.closest('table')
    const colIndex = cell.cellIndex
    const rowIndex = cell.parentNode.rowIndex
    if(colIndex === 15) return
    if(rowIndex<5 && colIndex<3 || rowIndex>4 && colIndex<1) return
    if(colIndex-(rowIndex>4 ? 1:3) === currentindex) return
    table.querySelectorAll("td").forEach(cell=>{
      cell.classList.remove("highlight")
    })
    table_x1.querySelectorAll("td").forEach(cell=>{
      cell.classList.remove("highlight")
    })
    Array.from(table.rows).forEach((row,i)=>{
      row.cells[colIndex-(i>4 ? 2:0)+(rowIndex>4 ? 2:0)].classList.add("highlight")
    })
    currentindex = colIndex-(rowIndex>4 ? 1:3)
    plot_hough(s,sigma,1,a_asym[currentindex])
  })
}

let p_r_s
async function loadLegendre(){
  const response = await fetch('nalegendre.bin')
  const buffer = await response.arrayBuffer()
  p_r_s = new Float64Array(buffer)
}

const canvas = document.getElementById("canvas")
const margin = {bottom:10,top:10,left:20,right:50}
const plot_height = canvas.height - margin.left - margin.top
const plot_width = canvas.width - margin.left - margin.right
const scaleX = (x) => x*plot_width/180
const scaleY = (y) => -y*plot_height/6
const ctx = canvas.getContext('2d')
ctx.translate(plot_width/2+margin.left,plot_height/2+margin.top)
ctx.font = '20px sans-serif'
ctx.textAlign = 'center'
ctx.textBaseline = 'middle'
ctx.save()

function plot_hough(s,sigma,flag_asymmetric,coef){
  const nlat = 91
  const m_max = 90
  let hough = new Float64Array(nlat)
  for (let i=flag_asymmetric;i<coef.length;i++){//?
    const start = (s*m_max+2*i-(flag_asymmetric===0 ? 0:1))*nlat
    const slice = p_r_s.subarray(start,start+nlat)
    for (let j=0;j<nlat;j++){
      hough[j]+= coef[flag_asymmetric===0 ? i:i-1]*slice[j]
    }
    //console.log(i,coef[i-1],slice)
  }
  const hough_t = Array.from(hough.subarray(1).map(e => e*(flag_asymmetric === 0 ? 1:-1))).reverse().concat(Array.from(hough))
  const hough_u = []
  const hough_v = []
  const mu = Array.from(Array(nlat*2-1),(_,i)=>Math.sin(Math.PI*(i-90)/180))

  let dtdm = (hough_t[1]-0)/(mu[1]-mu[0])
  let factor = Math.sqrt(1-mu[1]**2)/(sigma**2-mu[1]**2)
  hough_u[0] = factor*(s/(1-mu[1]**2)*hough_t[1] - mu[1]/sigma*dtdm)
  hough_v[0] = factor*(s*mu[1]/(sigma*(1-mu[1]**2))*hough_t[1] - dtdm)
  for (let i=1;i<hough_t.length-1;i++){
    const dtdm = (hough_t[i+1]-hough_t[i-1])/(mu[i+1]-mu[i-1])
    const factor = Math.sqrt(1-mu[i]**2)/(sigma**2-mu[i]**2)
    hough_u[i] = factor*(s/(1-mu[i]**2)*hough_t[i] - mu[i]/sigma*dtdm)
    hough_v[i] = factor*(s*mu[i]/(sigma*(1-mu[i]**2))*hough_t[i] - dtdm)
  }
  dtdm = (0-hough_t[hough_t.length-2])/(mu[hough_t.length-1]-mu[hough_t.length-2])
  factor = Math.sqrt(1-mu[hough_t.length-2]**2)/(sigma**2-mu[hough_t.length-2]**2)
  hough_u[hough_t.length-1] = factor*(s/(1-mu[hough_t.length-2]**2)*hough_t[hough_t.length-2] - mu[hough_t.length-2]/sigma*dtdm)
  hough_v[hough_t.length-1] = factor*(s*mu[hough_t.length-2]/(sigma*(1-mu[hough_t.length-2]**2))*hough_t[hough_t.length-2] - dtdm)

  hough_u[60] = (hough_u[59]+hough_u[61])/2
  hough_v[60] = (hough_v[59]+hough_v[61])/2
  hough_u[120] = (hough_u[119]+hough_u[121])/2
  hough_v[120] = (hough_v[119]+hough_v[121])/2
  document.getElementById('Hough_T').innerText = hough_t.join(', ')
  document.getElementById('Hough_U').innerText = hough_u.join(', ')
  document.getElementById('Hough_V').innerText = hough_v.join(', ')
  const norm_factor = parseFloat((Math.max(...hough_u.map(n=>Math.abs(n)))/Math.max(...hough_t.map(n=>Math.abs(n)))).toExponential(0))

  ctx.clearRect(-plot_width/2-margin.left,-plot_height/2-margin.top,canvas.width,canvas.height)
  ctx.restore()
  ctx.strokeStyle = 'black'
  ctx.lineWidth = 2
  ctx.strokeRect(scaleX(-90),scaleY(-3),scaleX(180),scaleY(6))
  for (let i=-90;i<=90;i+=10){
    ctx.beginPath()
    ctx.moveTo(scaleX(i),scaleY(-3))
    ctx.lineTo(scaleX(i),scaleY(-3.05))
    ctx.stroke()
  }
  for (let i=-3;i<=3;i+=0.2){
    ctx.beginPath()
    ctx.moveTo(scaleX(-90),scaleY(i))
    ctx.lineTo(scaleX(-90.5),scaleY(i))
    ctx.stroke()
  }
  ctx.lineWidth = 0.5
  ctx.setLineDash([3,3])
  for (let i=-60;i<90;i+=30){
    ctx.beginPath()
    ctx.moveTo(scaleX(i),scaleY(-3))
    ctx.lineTo(scaleX(i),scaleY(3))
    ctx.stroke()
    ctx.fillText(String(i).padStart(3,' ')+'°',scaleX(i),scaleY(-3.3))
  }
  for (let i=-2;i<3;i+=1){
    ctx.beginPath()
    ctx.moveTo(scaleX(-90),scaleY(i))
    ctx.lineTo(scaleX(90),scaleY(i))
    ctx.stroke()
    ctx.textAlign = "right"
    ctx.fillText(String(i),scaleX(-91),scaleY(i))
    ctx.textAlign = "left"
    ctx.fillStyle = 'red'
    ctx.fillText(String(i*norm_factor),scaleX(91),scaleY(i))
    ctx.fillStyle = 'black'
  }

  ctx.lineWidth = 3
  ctx.setLineDash([])


  if(document.getElementById('hough_t').checked){
    ctx.beginPath()
    ctx.moveTo(scaleX(-90),scaleY(hough_t[0]))
    for (let i=1;i<hough_t.length;i++){
      ctx.lineTo(scaleX(i-90),scaleY(hough_t[i]))
    }
    ctx.stroke()
  }
  if(document.getElementById('hough_u').checked){
    ctx.strokeStyle = 'red'
    ctx.beginPath()
    ctx.moveTo(scaleX(-90),scaleY(hough_u[0]))
    for (let i=1;i<hough_u.length;i++){
      ctx.lineTo(scaleX(i-90),scaleY(hough_u[i]/norm_factor))
    }
    ctx.stroke()
  }
  if(document.getElementById('hough_v').checked){
    ctx.strokeStyle = 'blue'
    ctx.beginPath()
    ctx.moveTo(scaleX(-90),scaleY(hough_v[0]))
    for (let i=1;i<hough_v.length;i++){
      ctx.lineTo(scaleX(i-90),scaleY(hough_v[i]/norm_factor))
    }
    ctx.stroke()
  }

}
window.addEventListener('load',async () => {
  await loadLegendre();
})

