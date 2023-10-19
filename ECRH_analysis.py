
from func_load import *

simloc = getSimulation('/storage/space2/phrmsf/ECRH/ECRH_zeroFE')#JT60U_5_2')
times = read_pkl('times')
d0 = sdfread(0)
klim = 0.5*2*const.PI/getdxyz(d0)
wlim = 0.5*2*const.PI/getdt(times)
wnorm = getCyclotronFreq(d0,'Protons')
fcp = getCyclotronFreq(d0,'Protons')/(2*const.PI)
tcp = 1/fcp
vA = getAlfvenVel(d0)
knorm = wnorm/vA
quantity = 'Magnetic_Field_Bz'
theta = 86.3*const.PI/180

## particle energies
# use rel_elec.py in energy/

# FT2d & power
wmax = 35
kmax = 35
FT2d = read_pkl('FT_2d_'+quantity)
nw,nk = FT2d.shape
FT2d = FT2d[:int(nw*wmax*wnorm/wlim),:int(nk*kmax*knorm/klim)]
nw,nk = FT2d.shape
#log10_power,omegas=powerspectrum(FT2d,wnorm,wmax,kmax,0,wmax,0,kmax)
#plt.plot(omegas/wnorm,log10_power,color='r',linestyle='--')
#print(len(omegas),len(log10_power))
#plt.show()
#sys.exit()
dk = klim/nk
dw = wlim/nw
print(dk/knorm,dw/wnorm)

## plotting 
#fig,axs = plt.subplots(ncols=2,figsize=(10,5),gridspec_kw={'width_ratios':[3,1]})#,sharey=True)
#fig.subplots_adjust(wspace=0.)
#ax=axs.ravel()
## ft2d
#ax[0].imshow(np.log10(FT2d[1:,1:]),extent=[0,kmax,0,wmax],**kwargs,cmap='magma') # clim=(-2,6)
#ax[0].set_xlim(0,kmax)
#ax[0].set_ylim(0,wmax)
## cold plasma disp
#omegas = wnorm*np.linspace(0,wmax,1000)
#k = omegas/(vA*np.sqrt(1+np.cos(theta)**2)) # Sumida 2023 eq. (3)
#ax[0].plot(k/knorm,omegas/wnorm,color='k',linestyle='--')
#ax[0].set_xlabel(r'$kv_A/\Omega_p$',fontsize=20) # \lambda_{De}
#ax[0].set_ylabel(r'$\omega/\Omega_p$',fontsize=20)
#ax[0].set_xticklabels([0,5,10,15,20,25,30,35,40]) #[0,0.005,0.01,0.015,0.02,0.025,0.03,0.035]
#ax[0].tick_params(axis='both',direction='out',top=False,right=False,left=True,bottom=True)
#log10_power,omegas=powerspectrum(FT2d,wnorm,[wmax,kmax],[0,wmax,0,kmax])
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
#ax[1].set_ylabel(r'$\omega/\Omega_p$',fontsize=18)
#fig.savefig('/storage/space2/phrmsf/ECRH/paper/FT2d_Bz_PSD.png',bbox_inches='tight')
#fig.savefig('/storage/space2/phrmsf/ECRH/paper/FT2d_Bz_PSD.eps',bbox_inches='tight')
#plt.show()

## bicoherence
# finer granularity (higher nffts)
#fieldmatrix = load_batch_fieldmatrix([],quantity)
#karea = klim/knorm #80
#warea = wlim/wnorm # 60
#nfft = len(times)//5
#noverlap = nfft//2
#bispec=True #boolean, True as default will calculate both bicoh AND bispec
#klabel = r'$v_A/\Omega_e$'
#fig, ax = getBicoh(karea,warea,fieldmatrix,getdt(times),times[-1],getGridlen(sdfread(0)),wnorm,knorm,nfft=nfft,\
#	noverlap=noverlap,window=True,bispectrum=bispec,klabel=klabel)


## spectrograms
#nfo = [10,2]
#try:
#	Sarr, freqs, time = read_pkl('Sxftarray_nfft_1645_noverlap_822') # 3291, 1645
#	Spower = read_pkl('Sxpower_nfft_1645_noverlap_822')
#except:
#	fieldmatrix = load_batch_fieldmatrix([],quantity)
#	Spower, Sarr, freqs, time = getSpectrogram(fieldmatrix,times,majspec='Protons',nfo=[10,2],minspec='FElectrons',plot=False)
#print(Sarr.shape)
#fig,ax = plt.subplots(nrows=3,figsize=(6,10))
#fig.subplots_adjust(hspace=0.1)
#nx,nf,nt = Sarr.shape
#upper=[nx//3,2*nx//3,nx-1]
#lower=[0,nx//3,2*nx//3]
#label=['(a)','(b)','(c)']
#i=0
#for l,u in zip(lower,upper): #[nx//3,2*nx//3,nx-1]: #range(0,nx,1000):
#	tSarr = np.mean(Sarr[l:u,:,:],axis=0)	
#	im = ax[i].imshow(tSarr,extent=[0,time[-1]/tcp,freqs[0],freqs[-1]/fcp],aspect='auto',interpolation='hanning',cmap='jet',clim=(0,1.7e-16)) # Sarr[x,:,:]
#	#plt.colorbar(im)
#	ax[i].set_ylim(0,100)
#	ax[i].annotate(label[i],xy=(0.9,0.1),xycoords='axes fraction',color='white',**tnrfont)
#	ax[i].set_ylabel(r'$f/f_{cp}$',**tnrfont)	
#	i+=1
#ax[-1].set_xlabel(r'$t\Omega_{p}/2\pi$',**tnrfont)
#fig.savefig('Sxft_meanoverX.png')
#plt.show()

## empirical growth rates
## change in du
#tplot, graddu = grad_energydens(simloc,normspecies='Protons',quant='Electrons',conv2=True)
#plt.plot(tplot,graddu)
#plt.show()
# FT1d growth rates
omegas, growthRatesMean, growthRatesSTD = map_k_growth(simloc,'Protons',0,30,tstart_frac=0.0,tend_frac=1.0,color='k')
thresh = growthRatesMean > 0
growthRatesMean = growthRatesMean[thresh]
omegas = omegas[thresh]
plt.plot(omegas/(2*const.PI*fcp),growthRatesMean/(2*const.PI*fcp),'-o')
plt.show()
