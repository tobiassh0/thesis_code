
from func_load import *

## load sim
simloc = getSimulation('/storage/space2/phrmsf/traceT/90deg')#lowres_D_He3/0_10_p_90')
d0 = sdfread(0)
times = read_pkl('times')

## set normalisation and limits
minspec = 'Alphas'
wnorm = getCyclotronFreq(d0,minspec)
vA = getAlfvenVel(d0)
kmax = 0.5*2*const.PI/getdxyz(d0)
wmax = 0.5*2*const.PI/getdt(times)
knorm = wnorm/vA
omegas,log10_power=power(wnorm,wklims=[wmax/wnorm,kmax/knorm],wkmax=[40,100],norm_omega=getOmegaLabel(minspec),quantity='Magnetic_Field_Bz',plot=False,read=False,dump=False,outp=True)
dw = 2*wmax/len(times) # account for +ve and -ve wmax covered by N_t files

#omegas = read_pkl('omegas_power')
power = 10**log10_power #read_pkl('log10_power')
thresh = omegas/wnorm < 25
plt.plot(omegas[thresh]/wnorm,power[thresh]/dw)
plt.ylabel('PSD',**tnrfont)
plt.yscale('log')
plt.xlabel(r'$\omega/$'+getOmegaLabel(minspec),**tnrfont)
plt.savefig('/home/space/phrmsf/Documents/thesis_code/power.png')
