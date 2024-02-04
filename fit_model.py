import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as opt
import pickle

import dsigma as ds
import parameters as pm

plt.style.use('standard')

data=np.loadtxt('fpi_data.csv',unpack=True,delimiter=',')
data_amendolia=np.loadtxt('fpi_amendolia.csv',unpack=True,delimiter=',')

nboot=10
Fpi=[]
g1=[]
Qs=[]
nsets=6
ndata=5
fit_data=False
if fit_data==True:
    for ii in range(nsets):
        qs=-data[0][ii*ndata]
        Qs.append(-qs)
        s=data[1][ii*ndata]**2
        t=-data[4][ii*ndata:(ii+1)*ndata]

        dL=data[5][ii*ndata:(ii+1)*ndata]
        dL_err=np.sqrt(data[6][ii*ndata:(ii+1)*ndata]**2+data[7][ii*ndata:(ii+1)*ndata]**2)
        dT=data[8][ii*ndata:(ii+1)*ndata]
        dT_err=data[9][ii*ndata:(ii+1)*ndata]**2
        dLT=data[10][ii*ndata:(ii+1)*ndata]
        dLT_err=data[11][ii*ndata:(ii+1)*ndata]**2
        dTT=data[12][ii*ndata:(ii+1)*ndata]
        dTT_err=data[13][ii*ndata:(ii+1)*ndata]**2

        boots=ds.bootstrap(dL,dL_err,nboot)
        boots=np.append(boots,[dL],axis=0)

        fpi_buf=[]
        g1_buf=[]
        # Fitting data
        for boot in boots:
            # Model has 2 free parameters, Fpi and g1.
            # They are passed as an array [Fpi, g1]
            fit=opt.curve_fit(lambda t_,fpi_,g1_: ds.dsLdt(s,t_,qs,[fpi_,g1_]),t[1:],boot[1:],sigma=boots.std(0)[1:])
            fpi_buf.append(fit[0][0])
            g1_buf.append(fit[0][1])
        Fpi.append(fpi_buf)
        g1.append(g1_buf)
    Fpi=np.array(Fpi).transpose()
    g1=np.array(g1).transpose()

    f=open('Fpi.dat','w')
    pickle.dump(Fpi,f,protocol=2)
    f.close()
    
    f=open('g1.dat','w')
    pickle.dump(g1,f,protocol=2)
    f.close()
    
else:
    f=open('Fpi.dat','r')
    Fpi=pickle.load(f)
    f.close()
    
    f=open('g1.dat','r')
    g1=pickle.load(f)
    f.close()

Qs=[]
for ii in range(nsets):
    qs=-data[0][ii*ndata]
    Qs.append(-qs)
    s=data[1][ii*ndata]**2
    t=-data[4][ii*ndata:(ii+1)*ndata]
    
    dL=data[5][ii*ndata:(ii+1)*ndata]
    dL_err=np.sqrt(data[6][ii*ndata:(ii+1)*ndata]**2+data[7][ii*ndata:(ii+1)*ndata]**2)
    dT=data[8][ii*ndata:(ii+1)*ndata]
    dT_err=data[9][ii*ndata:(ii+1)*ndata]**2
    dLT=data[10][ii*ndata:(ii+1)*ndata]
    dLT_err=data[11][ii*ndata:(ii+1)*ndata]**2
    dTT=data[12][ii*ndata:(ii+1)*ndata]
    dTT_err=data[13][ii*ndata:(ii+1)*ndata]**2
    
    boots=ds.bootstrap(dL,dL_err,nboot)
    boots=np.append(boots,[dL],axis=0)
    
    lam=[Fpi[-1,ii],g1[-1,ii]]
    print lam
    
    t_hd=np.linspace(min(t),max(t))
    
    # Predictions for cross sections
    fig,ax=plt.subplots(1,1)
    ax.plot(-t_hd,ds.dsLdt(s,t_hd,qs,lam))
    ax.errorbar(-t,dL,yerr=dL_err,marker='.',ls='')
    ax.errorbar(-t,boots.mean(0),yerr=boots.std(0),marker='x',ls='')
    ax.set_xlabel('$-t$ (GeV$^2$)')
    ax.set_ylabel(r'$d\sigma_\text{L}/dt$')
    fig.tight_layout()
    
    fig,ax=plt.subplots(1,1)
    ax.plot(-t_hd,ds.dsTdt(s,t_hd,qs,lam))
    ax.errorbar(-t,dT,yerr=dT_err,marker='.',ls='')
    ax.set_xlabel('$-t$ (GeV$^2$)')
    ax.set_ylabel(r'$d\sigma_\text{T}/dt$')
    fig.tight_layout()
    
    fig,ax=plt.subplots(1,1)
    ax.plot(-t_hd,ds.dsTTdt(s,t_hd,qs,lam))
    ax.errorbar(-t,dTT,yerr=dTT_err,marker='.',ls='')
    ax.set_xlabel('$-t$ (GeV$^2$)')
    ax.set_ylabel(r'$d\sigma_\text{TT}/dt$')
    fig.tight_layout()
    
    fig,ax=plt.subplots(1,1)
    ax.plot(-t_hd,ds.dsLTdt(s,t_hd,qs,lam))
    ax.errorbar(-t,dLT,yerr=dLT_err,marker='.',ls='')
    ax.set_xlabel('$-t$ (GeV$^2$)')
    ax.set_ylabel(r'$d\sigma_\text{LT}/dt$')
    fig.tight_layout()


# Plotting resulting extracted values for Fpi
Qs_hd=np.linspace(0,3)
fig,ax=plt.subplots(1,1)
fm2GeV=1/0.197 # GeV^-1
rPion=0.672*fm2GeV
ax.plot(Qs_hd,ds.fdip(Qs_hd,np.sqrt(6/rPion**2)))
ax.errorbar(Qs,Fpi[-1,:],Fpi.std(0),marker='.',ls='')
ax.errorbar(data_amendolia[0],data_amendolia[1],yerr=data_amendolia[2],marker='.',ls='',color='k',lw=1,capsize=2)

ax.set_xlabel('$Q^2$ (GeV$^2$)')
ax.set_ylabel('$F_\pi(Q^2)$')

ax.set_ylim([0,1])
ax.set_xlim([0,3])

fig.tight_layout()
plt.show()
