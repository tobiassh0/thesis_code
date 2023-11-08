
from func_load import *
import correlation.spatial_crosscor as sc


wca = 2.1 * (2 * abs(const.qe))/getMass('Alphas')
w = np.arange(0,40,0.01) * wca

simloc = getSimulation('/storage/space2/phrmsf/traceT/traceT_0_00')
#dphase,Rtdx = sc.getSpatialCorrelation(sdfread(0),plot=False)
vA = getAlfvenVel(sdfread(0))

## cold plasma disp frequencies
_,k2,_= coldplasmadispersion(sdfread(0),w,theta=89.)
k2 *= (vA/wca)
plt.xlim(0,100)
plt.scatter(k2,w/wca,facecolor='b',edgecolor='b')

## 1994 oblique (24)
theta = 89. * const.PI/180.
k = np.arange(0,1000,0.01) * (wca/vA)
kpara = k*np.cos(theta)
w2 = (k**2 + kpara**2 + (k*kpara)**2 * (vA/wca)**2 )**2 - (2*k*kpara)**2
tw2 = 0.5*(vA**2)*(k**2 + kpara**2 + (k*kpara)**2 * (vA/wca)**2 + np.sqrt(w2))
tw = np.sqrt(tw2)
plt.scatter(k*vA/wca,tw/wca,facecolor='r',edgecolor='r')
plt.show()
