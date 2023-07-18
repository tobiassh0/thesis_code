
from func_load import *


sim_loc = getSimulation('/storage/space2/phrmsf/cold_JET26148')
#fm = load_batch_fieldmatrix()
#dfm = fm - np.mean(fm[0:10,:]) # delta Bz
times = read_pkl('times')

#FT2d = get2dTransform(dfm)
#dumpfiles(FT2d,'FT_2d_Magnetic_Field_Bz')
FT2d = read_pkl('FT_2d_Magnetic_Field_Bz')
wcyc = getCyclotronFreq(sdfread(0),'Deuterons')
vA = getAlfvenVel(sdfread(0))
wnorm = wcyc
knorm = wcyc/vA#1/getDebyeLength(sdfread(0),'Electrons')
klim_prime = (0.5*2*const.PI/getdxyz(sdfread(0)))/knorm
wlim_prime = (0.5*2*const.PI/getdt(times))/wnorm
wmax = 26
kmax = 40 #0.05

fig,axs = plt.subplots(ncols=2,figsize=(10,5),gridspec_kw={'width_ratios':[3,1]})#,sharey=True)
fig.subplots_adjust(wspace=0.)
ax=axs.ravel()
ax[0].imshow(np.log10(FT2d[1:]),extent=[0,klim_prime,0,wlim_prime],**kwargs,clim=(-2,6),cmap='magma')
ax[0].set_xlim(0,kmax)
ax[0].set_ylim(0,wmax)
ax[0].set_xlabel(r'$kv_A/\Omega_D$',fontsize=20) # \lambda_{De}
ax[0].set_ylabel(r'$\omega/\Omega_D$',fontsize=20)
ax[0].set_xticklabels([0,5,10,15,20,25,30,35,40]) #[0,0.005,0.01,0.015,0.02,0.025,0.03,0.035]
ax[0].tick_params(axis='both',direction='out',top=False,right=False,left=True,bottom=True)
omegas, log10_power = power(klim_prime,wlim_prime,wmax,kmax,wnorm,norm_omega=r'$\Omega_D$',quantity='Magnetic_Field_Bz',plot=True,read=False,outp=True)
dw = (omegas[-1]-omegas[0])/len(omegas)
psd = (10**log10_power)/dw
for i in np.arange(0,30):
	ax[1].axhline(i,color='darkgray',linestyle='--')
ax[1].plot(psd,omegas/wnorm)
ax[1].set_xscale('log')
ax[1].set_xlim(10**(-2.5),10**(0.5))
ax[1].set_ylim(0,wmax)
ax[1].set_xlabel(r'PSD',**tnrfont)#fontsize=20)
ax[1].set_yticklabels([])
ax[1].tick_params(axis='both',direction='out',top=False,right=False,left=True,bottom=True)

#ax[1].set_ylabel(r'$\omega/\Omega_D$',fontsize=18)

fig.savefig('/storage/space2/phrmsf/paper/Cold_FT2d_PSD.png',bbox_inches='tight')
fig.savefig('/storage/space2/phrmsf/paper/Cold_FT2d_PSD.eps',bbox_inches='tight')
plt.show()


