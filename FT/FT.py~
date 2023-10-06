from func_load import *

simloc = getSimulation('/storage/space2/phrmsf/ECRH/ECRH_JT60U_5_2')
d0 = sdfread(0)
times = read_pkl('times')
klim = 0.5*2*const.PI/getdxyz(d0)
wlim = 0.5*2*const.PI/getdt(times)
wce = getCyclotronFreq(d0,'Electrons')
wcp = getCyclotronFreq(d0,'Protons')
vA = getAlfvenVel(d0)
knorm = wce/vA
tlim = times[-1]
tce = 2*const.PI/wce

### 1d
#FT1d = read_pkl('FT_1d_Magnetic_Field_Bz')
#fig,ax=plot1dTransform(FT1d,klim/knorm,tlim/tce,klabel=r'$v_A/\Omega_e$',wlabel=r'$\Omega_e$')
#fig.savefig('FT_1d_Magnetic_Field_Bz.png')

### 2d
#FT2d = read_pkl('FT_2d_Magnetic_Field_Bz')
#fig,ax=plot2dTransform(FT2d,klim/knorm,wlim/wce,klabel=r'$v_A/\Omega_e$',wlabel=r'$\Omega_e$')
#plotting(fig,ax,'FT_2d_Magnetic_Field_Bz')

### power
#power(klim*vA/wcp,wlim/wcp,35,100,wcp,norm_omega=r'$\Omega_p$',quantity='Magnetic_Field_Bz',plot=True,read=False,outp=True)
