
from func_load import *


simloc = getSimulation('/storage/space2/phrmsf/ECRH/ECRH_JT60U_5_2')
#FT2d = read_pkl('FT_2d_Magnetic_Field_Bz')
times = read_pkl('times')
d0 = sdfread(0)
kmax = 0.5*2*const.PI/getdxyz(d0)
wmax = 0.5*2*const.PI/getdt(times)
wnorm = getCyclotronFreq(d0,'Protons')
vA = getAlfvenVel(d0)
knorm = wnorm/vA

## load and plot
#plt.imshow(np.log10(FT2d[1:,1:]),**kwargs,extent=[0,kmax*vA/wnorm,0,wmax/wnorm])
#plt.show()
#
#omegas, log10power = power(kmax*vA/wcp,wmax/wcp,60,80,wnorm,norm_omega=r'$\Omega_p$',quantity='Magnetic_Field_Bz',plot=True,read=False)
#plt.plot(omegas/wnorm,log10power)
#plt.show()

## replot
simloc = getSimulation('/storage/space2/phrmsf/lowres_D_He3/0_10_p_90')
_,_=power(kmax/knorm,wmax/wnorm,40,100,wnorm,norm_omega=r'$\Omega_p$',quantity='Magnetic_Field_Bz',plot=True,read=True,outp=False)

omegas = read_pkl('omegas_power')
log10power = read_pkl('log10_power')
plt.xlim(0,40)
plt.plot(omegas/wnorm,log10power) ; plt.show()