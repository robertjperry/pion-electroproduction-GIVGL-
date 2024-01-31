import numpy as np
import matplotlib.pyplot as plt

import dsigma as ds
import parameters as pm

plt.style.use('standard')

s=1.95**2 # GeV^2
t=-0.2 # GeV^2
qs=-0.6 # GeV^2
#qs=0.14**2
p1,q,p2,k=ds.get_fv(s,t,qs)
print p1
print q
print p2
print k

print p1+q
print p2+k
print ds.dot(q-k,q-k)
print ds.dot(p1-p2,p1-p2)

d=q-k
print d[0]**2,np.dot(d[1:],d[1:])

#quit()

t=np.linspace(-0.08,-0.04)

lam=[ds.kF1p(qs),0.0]

fig,ax=plt.subplots(1,1)
ax.plot(-t,np.real(ds.dsLdt(s,t,qs,lam)))
ax.plot(-t,np.imag(ds.dsLdt(s,t,qs,lam)))
ax.set_xlabel('$-t$ (GeV$^2$)')
ax.set_ylabel(r'$d\sigma_\text{L}/dt$')
fig.tight_layout()

fig,ax=plt.subplots(1,1)
ax.plot(-t,np.real(ds.dsTdt(s,t,qs,lam)))
ax.plot(-t,np.imag(ds.dsTdt(s,t,qs,lam)))
ax.set_xlabel('$-t$ (GeV$^2$)')
ax.set_ylabel(r'$d\sigma_\text{T}/dt$')
fig.tight_layout()

fig,ax=plt.subplots(1,1)
ax.plot(-t,np.real(ds.dsTTdt(s,t,qs,lam)))
ax.plot(-t,np.imag(ds.dsTTdt(s,t,qs,lam)))
ax.set_xlabel('$-t$ (GeV$^2$)')
ax.set_ylabel(r'$d\sigma_\text{TT}/dt$')
fig.tight_layout()

fig,ax=plt.subplots(1,1)
ax.plot(-t,np.real(ds.dsLTdt(s,t,qs,lam)))
ax.plot(-t,np.imag(ds.dsLTdt(s,t,qs,lam)))
ax.set_xlabel('$-t$ (GeV$^2$)')
ax.set_ylabel(r'$d\sigma_\text{LT}/dt$')
fig.tight_layout()



plt.show()
