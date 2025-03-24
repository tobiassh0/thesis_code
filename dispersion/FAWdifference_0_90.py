
from func_load import *

B0 = 2.07
n0 = 2.4611176521760727e19
n90 = 1.7e19
species = ['Electrons','Deuterons','Tritons','Alphas']
charges = [const.qe*getChargeNum(spec) for spec in species]
masses = [getMass(spec) for spec in species]
wcf = [charges[i]*B0/masses[i] for i in range(len(masses))]
xi3 = 1.5e-4
xiarr0 = [1, 1-2*xi3, 0 , xi3]
xiarr90 = [1, 1-0.9-2*xi3, 0.9, xi3]

# VA - should be equivalent
vA0=B0/np.sqrt(const.mu0*np.sum([xiarr0[i]*masses[i]*n0 for i in range(len(masses))]))
vA90=B0/np.sqrt(const.mu0*np.sum([xiarr90[i]*masses[i]*n90 for i in range(len(masses))]))
print(vA0,vA90,vA0/vA90)

wpf0 = [np.sqrt((n0*xiarr0[i]*charges[i]**2)/(masses[i]*const.e0)) for i in range(len(masses))]
wpf90 = [np.sqrt((n90*xiarr90[i]*charges[i]**2)/(masses[i]*const.e0)) for i in range(len(masses))]
_theta = 90*const.PI/180
omegas = np.linspace(0,40,1000)*wcf[-1]

# 0%
k1,_k2,k3 = coldplasmadispersion_analytical(omegas,wpf0,wcf,theta=_theta)
plt.plot(np.real(_k2)*vA0/wcf[-1],omegas/wcf[-1],color='b',label=r'$0\%$')
print(np.real(_k2))

# 90%
k1,k2,k3 = coldplasmadispersion_analytical(omegas,wpf90,wcf,theta=_theta)
plt.plot(np.real(k2)*vA90/wcf[-1],omegas/wcf[-1],color='r',label=r'$90\%$')

# frequency difference
index0 = np.argmin(np.abs(15-_k2*vA0/wcf[-1]))
index90 = np.argmin(np.abs(15-k2*vA90/wcf[-1]))
plt.scatter(_k2[index0]*vA0/wcf[-1],omegas[index0]/wcf[-1],facecolor='b')
plt.scatter(k2[index90]*vA90/wcf[-1],omegas[index90]/wcf[-1],facecolor='r')
print(_k2[index0]*vA0/wcf[-1],omegas[index0]/wcf[-1])

plt.xlim(0,35)#10,20)
plt.ylim(0,25)#15,35)
plt.axvline(15,linestyle='--',color='k')
plt.legend(loc='best')
plt.ylabel(r'$\omega/\Omega_\alpha$',**tnrfont)
plt.xlabel(r'$k v_A/\Omega_\alpha$',**tnrfont)

plt.savefig('FAW_disp_frequency_diff_theta_{:.1f}.png'.format(_theta*180/const.PI),bbox_inches='tight')
plt.show()
