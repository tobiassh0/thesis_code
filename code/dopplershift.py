from func_load import *
#from list_new import *

### General process of doppler shift is a linear trend if the frequency (y-axis) due to the wavenumber (x-axis)
## un-commenting will reveal this graphically with a simple example  
#x = np.arange(0,100,0.01)
#for i in range(0,30,2):
#	y = np.ones(len(x))*i
#	yds = y - x
#	plt.plot(x,yds,color='k',alpha=1-i/30)
#plt.show()

# load sim & FT2d 
sim_loc = getSimulation('/storage/space2/phrmsf/lowres_D_He3/0_25_p_90')
FT2d = read_pkl('FT_2d_Magnetic_Field_Bz')
d0 = sdfread(0)
times = read_pkl('times')
vA = getAlfvenVel(d0)
print(vA)
Ep = 14.68*1e6*const.qe
vp = (2*Ep/getMass('Protons'))**0.5 # proton velocity
klim = 0.5*2*const.PI/getdxyz(d0)
wlim = 0.5*2*const.PI/getdt(times)
wcp = getCyclotronFreq(d0,'Protons')
wcD = getCyclotronFreq(d0,'Deuterons')
wnorm = wcp
knorm = wnorm/vA
dk = 5*2*const.PI/getGridlen(d0)
dw = 0.5*2*const.PI/times[-1]

## calc gradients in image
# cut FT2d into size needed
(nk,nw) = FT2d.shape
kmax = 50*knorm ; wmax = 10*wnorm
FT2d = np.log10(FT2d[:int(nw*wmax/wlim),:int(nk*kmax/klim)])
(nw,nk) = FT2d.shape
# threshold FT2d
thresh = 1.8
for i in range(nw):
	for j in range(nk):
		if FT2d[i,j] < thresh:
			FT2d[i,j] = 0
		else:
			continue
# Scharr gradient map
_,kGangle = Kernel(FT2d,kernel='sobel') # scharr or sobel
# gradients as angles
kGangle = kGangle[1:-1,1:-1]# remove abberations around edge
kGangle = kGangle.flatten()
# convert to all negative angles (easier to calc real gradient)
for i in range(len(kGangle)):
	if kGangle[i] > 0:
		kGangle[i]-=const.PI #rad
# remove inf values
dwdk = np.nan_to_num(np.tan(kGangle),posinf=np.nan,neginf=np.nan) # remove inf and -inf values
dwdk = dwdk[~np.isnan(dwdk)]
# remove zero values (dont want to plot them in hist)
dwdk = dwdk[dwdk!=0]
dw_dk = dwdk * (dw/dk)/vA # normalise
thresh = (np.abs(dw_dk) < 10)
dw_dk = dw_dk[thresh]
print('Scharr kernel mean :: ',np.mean(dw_dk))
print('Scharr kernel medi :: ',np.median(dw_dk))
# plot hist
counts,bins,_=plt.hist(dw_dk,bins=10000,range=(-1,1),density=True) # np.log10
print('Scharr kernel max :: ', bins[np.argmax(counts)])
plt.ylabel('Normalised count',**tnrfont)
plt.xlabel(r'$d\omega/dk$'+'  '+r'$[v_A]$',**tnrfont)
#plt.savefig('dw_dk_Sobel_grad.png')
plt.show()
sys.exit()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# calculate velocities and est. dopp shift
E0 = 14.57*1e6*const.qe # MeV to Joules
vp = np.sqrt(2*E0/getMass('Protons'))
vperp = 0.9*vA
vpara = (vp**2 - vperp**2)**.5
print(vA,vp,vA-vp)
ds = vA*((vA-vp)/const.c)
dds = vA*((vA-vpara)/const.c)
kx = np.linspace(0,100,1000)*knorm
wmax = 20
#for i in range(0,wmax,1):
#	w = wcp*np.ones(len(kx))*i
#	ww = w + ds*kx
#	plt.plot(kx/knorm,ww/wcp,'k--')
##	ww = w + dds*kx
##	plt.plot(kx/knorm,ww/wcp,'b--')
#plt.ylabel(r'$\omega/\Omega_D$',fontsize=18)
#plt.xlabel(r'$k\lambda_{De}$',fontsize=18)
#plt.xlim(0,0.06)
#plt.ylim(0,wmax+5)
#plt.show()


vds = (vp-vA)*vA/const.c
(NW, NK) = FT2d.shape
print(NW,NK)
v = 0
dv = vds/100
## shift FT2d data and take power spectra
#tFT2d = np.zeros(FT2d.shape)
#for i in range(tFT2d.shape[1]):
#	ki = (i/NK)*klim # real valued
#	for j in range(FT2d.shape[0]):
#		wj = (j/NW)*wlim # real valued
#		wds = wj + ki*vds
#		jj = (wds/wlim)*NW ## undo doppler shift
#		try:
#			if tFT2d[int(jj),i] != 0: # already contains data
#				None
#			else:
#				tFT2d[int(jj),i] = FT2d[j,i] # shifted data
#		except:
#			None
## plot shifted FT2d
#fig, ax = plt.subplots()
#ax[1].imshow(np.log10(tFT2d),**kwargs,extent=[0,klim/knorm,0,wlim/wcp],cmap='magma',clim=(-2,6))
#ax[1].set_yticklabels([])
#ax[1].set_xlabel(r'$k\lambda_{De}$',fontsize=18)
#ax[0].set_ylim(0,35) ; ax[1].set_ylim(0,35)
#ax[0].set_xlim(0,0.06) ; ax[1].set_xlim(0,0.06)
#fig.savefig('/storage/space2/phrmsf/dump/dopplershift_FT2d.png')










