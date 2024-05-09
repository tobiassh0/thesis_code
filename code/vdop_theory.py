
from func_load import *

# setup & consts
theta_B = 89*const.PI/180
B0 = 3.7
n0 = 5e19
mp = const.me_to_mp*const.me
Ep = 14.68*1e6*const.qe
up = np.sqrt(2*Ep/mp)
#print(u0)

xiHe3 = np.array([0,0.05,0.1,0.15,0.22,0.25,0.34,0.38,0.45]) # np.linspace(0,0.45,100)
masses = [const.me_to_mD*const.me,const.me_to_He3*const.me,const.me_to_mp*const.me]
charges = [1,2,1]
xip = 1e-3
xiD = (1/charges[0])*(1-xiHe3*charges[1]-xip*charges[2])
rho = n0*(xiD*masses[0] + xiHe3*masses[1] + xip*masses[2])
vA = B0/np.sqrt(const.mu0*rho)
pitch = np.arcsin(0.9*vA/up)

vdop = (up/vA)*(np.cos(theta_B)*np.cos(pitch))
print(vdop)
plt.scatter(xiHe3,vdop)
plt.show()
