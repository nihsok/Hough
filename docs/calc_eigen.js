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
  for (let i=0;i<5;i++){
    let j = i==0 ? 2+i : i //rowspan for the first row
    row[i].cells[j].innerText = M[i*2].toFixed(3)
    row[i].cells[j+1].innerText = L[i*2].toFixed(3)
    if(i<4) row[i+1].cells[i].innerText = L[i*2].toFixed(3)
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
  let a_sym = []
  let n_p = s-2
  let n_n = s
  for (let i=1;i<=15;i++){
    x1[0].cells[i+2].innerHTML = 'Θ<sub>'+String(lambda[len_sym-i]>0 ? n_p+=2 : n_n-=2)+'</sub><sup>'+String(s)+','+String(sigma)+'</sup>'
    x1[1].cells[i+2].innerHTML = (4*Math.PI*7000/Math.sqrt(8000/(lambda[len_sym-i]*ev2ed)-1)).toFixed(1)
    x1[2].cells[i+2].innerText =(lambda[len_sym-i]*ev2ed).toFixed(1)
    x1[3].cells[i+2].innerText = lambda[len_sym-i].toFixed(4)
    if (lambda[len_sym-i]<0){
      x1[2].cells[i+2].innerHTML = '<i>'+x1[2].cells[i+2].innerText+'</i>'
      x1[3].cells[i+2].innerHTML = '<i>'+x1[3].cells[i+2].innerText+'</i>'
    }
    x1[4].cells[i+2].innerText = vector[len_sym-i].vector._data[0].toFixed(3)
    x1[5].cells[i].innerText = vector[len_sym-i].vector._data[1].toFixed(3)
    x1[6].cells[i].innerText = vector[len_sym-i].vector._data[2].toFixed(3)
    x1[7].cells[i].innerText = vector[len_sym-i].vector._data[3].toFixed(3)
    x1[8].cells[i].innerText = vector[len_sym-i].vector._data[4].toFixed(3)
    x1[9].cells[i].innerText = '⋮'

    a_sym[i-1] = vector[len_sym-i].vector._data
  }
  const table = document.getElementById('X1')
  let currentindex = null
  table.addEventListener("mouseover",function(e){
    const cell = e.target.closest("td")
    if(!cell) return
    const colIndex = cell.cellIndex
    const rowIndex = cell.parentNode.rowIndex
    if(colIndex-(rowIndex>4 ? 1:3) === currentindex) return
    if(rowIndex<5 && colIndex<3 || rowIndex>4 && colIndex<1) return
    if(currentindex !== null){
      Array.from(table.rows).forEach((row,i)=>{
        row.cells[currentindex+(i>4 ? 1:3)].classList.remove("highlight")
      })
    }
    currentindex = colIndex-(rowIndex>4 ? 1:3)
    Array.from(table.rows).forEach((row,i)=>{
      row.cells[colIndex-(i>4 ? 2:0)+(rowIndex>4 ? 2:0)].classList.add("highlight")
    })
    plot_hough(s,0,a_sym[currentindex])
  })
  //antisymmetric mode
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
const ctx = canvas.getContext('2d')
ctx.translate(plot_width/2+margin.left,plot_height/2+margin.top)
ctx.font = '20px sans-serif'
ctx.textAlign = 'center'
ctx.textBaseline = 'middle'
ctx.save()
const scaleX = (x) => x*plot_width/180
const scaleY = (y) => -y*plot_height/6

function plot_hough(s,flag_asymmetric,coef){
  //add title

  const nlat = 91
  const m_max = 90
  let hough = new Float64Array(nlat)
  for (let i=flag_asymmetric;i<n/2;i++){//n/2;i++){
    const start = s*m_max*nlat+2*i*nlat
    const slice = p_r_s.subarray(start,start+nlat)
    for (let j=0;j<nlat;j++){
      hough[j]+= coef[i]*slice[j]
    }
  }
  const hough_t = Array.from(hough.subarray(1)).reverse().concat(Array.from(hough))
  const hough_u = []
  const hough_v = []
  console.log(hough_t)

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

  ctx.restore()
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

