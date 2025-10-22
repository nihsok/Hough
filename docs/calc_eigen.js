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
    plot_hough(s,0,a_sym[currentindex])
  })
  table_x1.addEventListener('mouseleave', e=>{
    currentindex = null
    e.target.closest('table').querySelectorAll("td").forEach(cell=>{
      cell.classList.remove("highlight")
    })
  })


  //antisymmetric mode
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
    plot_hough(s,1,a_asym[currentindex])
  })
}

let p_r_s
async function loadLegendre(){
  const response = await fetch('nalegendre.bin')
  const buffer = await response.arrayBuffer()
  p_r_s = new Float64Array(buffer)
}

const canvas = document.getElementById("canvas")
const margin = {bottom:10,top:10,left:20,right:10}
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

function plot_hough(s,flag_asymmetric,coef){
  const nlat = 91
  const m_max = 90
  let hough = new Float64Array(nlat)
  for (let i=flag_asymmetric;i<coef.length;i++){//?
    const start = (s*m_max+2*i-(flag_asymmetric===0 ? 0:1))*nlat
    const slice = p_r_s.subarray(start,start+nlat)
    for (let j=0;j<nlat;j++){
      hough[j]+= coef[flag_asymmetric===0 ? i: i-1]*slice[j]
    }
    console.log(i,coef[i-1],slice)
  }
  const hough_t = Array.from(hough.subarray(1).map(e => e*(flag_asymmetric === 0 ? 1:-1))).reverse().concat(Array.from(hough))
  const hough_u = []
  const hough_v = []
  document.getElementById('Hough_T').innerText = hough_t.join(', ')

  ctx.clearRect(-plot_width/2-margin.left,-plot_height/2-margin.top,canvas.width,canvas.height)
  ctx.restore()
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
    ctx.fillText(String(i),scaleX(-94),scaleY(i))
  }

  ctx.lineWidth = 3
  ctx.setLineDash([])
  ctx.beginPath()
  ctx.moveTo(scaleX(-90),scaleY(hough_t[0]))
  for (let i=1;i<hough_t.length;i++){
    ctx.lineTo(scaleX(i-90),scaleY(hough_t[i]))
  }
  ctx.stroke()

}
window.addEventListener('load',async () => {
  await loadLegendre();
})

