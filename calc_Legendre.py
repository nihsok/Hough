#calculate normalized associated Legendre polynomials
from scipy.special import assoc_legendre_p_all
import numpy as np

m_max = 90 #order: m=m=s 50?
l_max= 2*m_max #degree: l=n=r

lp = assoc_legendre_p_all(l_max,m_max,np.sin(np.linspace(0,np.pi/2,91)),norm=True)
#print(lp[0,:,:,:].shape) #lp(0:0,0:n,-m:m,0:90)
with open("docs/nalegendre.bin",'wb') as f:
  for m in range(m_max):
    for l in range(m,m+m_max):
      f.write(lp[0,l,m,:])