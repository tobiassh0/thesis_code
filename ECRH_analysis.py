
from func_load import *

simloc = getSimulation('/storage/space2/phrmsf/ECRH/ECRH_JT60U_5_2')
times = read_pkl('times')
d0 = sdfread(0)
klim = 0.5*2*const.PI/getdxyz(d0)
wlim = 0.5*2*const.PI/getdt(times)
wnorm = getCyclotronFreq(d0,'Protons')
vA = getAlfvenVel(d0)
knorm = wnorm/vA
quantity = 'Magnetic_Field_Bz'

## FT2d & p#ower
#wmax = 35
#kmax = 100
#FT2d = read_pkl('FT_2d_'+quantity)
#nw,nk = FT2d.shape
#FT2d = FT2d[:int(nw*wmax*wnorm/wlim),:int(nk*kmax*knorm/klim)]
#nw,nk = FT2d.shape
##log10_power,omegas=powerspectrum(FT2d,wnorm,wmax,kmax,0,wmax,0,kmax)
##plt.plot(omegas/wnorm,log10_power,color='r',linestyle='--')
##print(len(omegas),len(log10_power))
##plt.show()
##sys.exit()
#dk = klim/nk
#dw = wlim/nw
#print(dk/knorm,dw/wnorm)

### plotting 
#fig,axs = plt.subplots(ncols=2,figsize=(10,5),gridspec_kw={'width_ratios':[3,1]})#,sharey=True)
#fig.subplots_adjust(wspace=0.)
#ax=axs.ravel()
#ax[0].imshow(np.log10(FT2d[1:,1:]),extent=[0,kmax,0,wmax],**kwargs,cmap='magma') # clim=(-2,6)
#ax[0].set_xlim(0,kmax)
#ax[0].set_ylim(0,wmax)
#ax[0].set_xlabel(r'$kv_A/\Omega_p$',fontsize=20) # \lambda_{De}
#ax[0].set_ylabel(r'$\omega/\Omega_p$',fontsize=20)
#ax[0].set_xticklabels([0,5,10,15,20,25,30,35,40]) #[0,0.005,0.01,0.015,0.02,0.025,0.03,0.035]
#ax[0].tick_params(axis='both',direction='out',top=False,right=False,left=True,bottom=True)
#log10_power,omegas=powerspectrum(FT2d,wnorm,wmax,kmax,0,wmax,0,kmax)
#dw = (omegas[-1]-omegas[0])/len(omegas)
#psd = (10**log10_power)/dw
#for i in np.arange(0,wmax):
#	ax[1].axhline(i,color='darkgray',linestyle='--')
#ax[1].plot(psd,omegas/wnorm)
#ax[1].set_xscale('log')
##ax[1].set_xlim(10**(-2.5),10**(0.5))
#ax[1].set_ylim(0,wmax)
#ax[1].set_xlabel(r'PSD',**tnrfont)#fontsize=20)
#ax[1].set_yticklabels([])
#ax[1].tick_params(axis='both',direction='out',top=False,right=False,left=True,bottom=True)
##ax[1].set_ylabel(r'$\omega/\Omega_D$',fontsize=18)
#fig.savefig('/storage/space2/phrmsf/ECRH/paper/FT2d_Ex_PSD.png',bbox_inches='tight')
#fig.savefig('/storage/space2/phrmsf/ECRH/paper/FT2d_Ex_PSD.eps',bbox_inches='tight')
#plt.show()

## bicoherence
# finer granularity (higher nffts)
fieldmatrix = load_batch_fieldmatrix([],quantity)
karea = 90
warea = 60
nfft = len(times)//5
noverlap = nfft//2
bispec=True #boolean, True as default will calculate both bicoh AND bispec
klabel = r'$v_A/\Omega_p$'
fig, ax = getBicoh(karea,warea,fieldmatrix,getdt(times),times[-1],getGridlen(sdfread(0)),wnorm,knorm,nfft=nfft,\
	noverlap=noverlap,window=True,bispectrum=bispec,klabel=klabel)

## spectrograms


## linear MCI growth rates
