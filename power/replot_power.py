
from func_load import *


simloc = getSimulation('/storage/space2/phrmsf/ECRH/ECRH_JT60U_5_2')
FT2d = read_pkl('FT_2d_Magnetic_Field_Bz')
times = read_pkl('times')
d0 = sdfread(0)
kmax = 0.5*2*const.PI/getdxyz(d0)
wmax = 0.5*2*const.PI/getdt(times)
wnorm = getCyclotronFreq(d0,'Protons')
vA = getAlfvenVel(d0)
plt.imshow(np.log10(FT2d[1:,1:]),**kwargs,extent=[0,kmax*vA/wnorm,0,wmax/wnorm])
plt.show()

omegas, log10power = power(kmax*vA/wcp,wmax/wcp,60,80,wnorm,norm_omega=r'$\Omega_p$',quantity='Magnetic_Field_Bz',plot=True,read=False)
plt.plot(omegas/wnorm,log10power)
plt.show()
