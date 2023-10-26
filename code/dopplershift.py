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
sim_loc = getSimulation('/storage/space2/phrmsf/lowres_D_He3/0_38_p_90')
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
(nw,nk) = (FT2d.shape)
dk = 2*const.PI/getGridlen(d0) # no factor half (see e-notes 24/10/23)
dw = 2*const.PI/times[-1]
print(dw,dk)
# figure setup
fig,ax=plt.subplots(ncols=3,figsize=(6*3,4))

## calc gradients in image
# cut FT2d into size needed
(nw,nk) = FT2d.shape
kmax = 20*knorm ; wmax = 10*wnorm
kmin = 0*knorm; wmin = 0#*wnorm
FT2d = np.log10(FT2d[:int(nw*wmax/wlim),int(nk*kmin/klim):int(nk*kmax/klim)])
(nw,nk) = FT2d.shape

# threshold FT2d
thresh = 1.8
tFT2d = FT2d.copy()
#tFT2d[FT2d > 1.6] = 0
tFT2d[FT2d < thresh] = 0

#for i in range(nw):
#	for j in range(nk):
#		
#		if FT2d[i,j] < thresh:
#			FT2d[i,j] = 0
#		else:
#			continue

# Kernel gradient map
kernel = 'custom'
_,kGangle = Kernel(tFT2d,kernel=kernel) # scharr, sobel or custom
# gradients as angles
kGangle = kGangle[1:-1,1:-1]# remove abberations around edge
im = ax[0].imshow(kGangle*180/const.PI,**kwargs,extent=[0,kmax/knorm,0,wmax/wnorm],cmap='Accent')
plt.colorbar(im)
print((dw/dk)/vA)
# remove zero values (dont want to plot them in hist)
dwdk = np.nan_to_num(np.tan(kGangle),posinf=np.nan,neginf=np.nan) # remove inf and -inf values
dw_dk = dwdk * (dw/dk)/vA # normalise
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
thresh = (np.abs(dw_dk) < 2.0) & (np.abs(dw_dk) > 0.001)
dw_dk = dw_dk[thresh]
print(kernel+' kernel mean :: ',np.mean(dw_dk))
print(kernel+' kernel medi :: ',np.median(dw_dk))
# plot hist
counts,bins,_=ax[1].hist(dw_dk,bins=1000,density=True,range=(-1,1)) # np.log10
dsv = bins[np.argmax(counts)] # doppler shift velocity in units of vA
print(kernel+' kernel max :: ', dsv)
ax[1].set_ylabel('Normalised count',**tnrfont)
ax[1].set_xlabel(r'$d\omega/dk$'+'  '+r'$[v_A]$',**tnrfont)
# plot kde
kde = stats.gaussian_kde(dw_dk)
xx = np.linspace(-1,1,1000)
ax[1].plot(xx,kde(xx),color='r')
#fig.savefig('dw_dk_'+kernel+'_grad.png',bbox_inches='tight')
#plt.show()

# plot FT2d
#fig,ax=plt.subplots(figsize=(8,6))
ax[2].imshow((tFT2d[1:,1:]),**kwargs,extent=[0,kmax/knorm,0,wmax/wnorm])
kx = np.linspace(0,20,100)*knorm
theta = 86.3 # deg
kperp = np.sin(theta*const.PI/180)
kpara = np.cos(theta*const.PI/180)
uperp = 0.9; upara = 6.076 
dsth = -(kperp*uperp + kpara*upara)
print(dsv,(dsth+1))
# doppler shifted line 
for i in range(0,int(wmax/wnorm),1):
	w = wnorm*np.ones(len(kx))*i
	# empirical
	ww = w + (dsv*vA)*kx
	ax[2].plot(kx/knorm,ww/wnorm,color='white',linestyle='--')
#	# theory
#	tww = w + ((dsth+1)*vA)*kx
#	ax.plot(kx/knorm,tww/wnorm,color='white',linestyle='-.')
ax[2].set_xlim(0,20)
ax[2].set_ylim(0,20)
ax[2].set_ylabel(r'$\omega/\Omega_p$',**tnrfont)
ax[2].set_xlabel(r'$kv_A/\Omega_p$',**tnrfont)
ax[2].plot([0,10],[0,10],color='white',linestyle='--') # vA line

#fig.savefig('FT_2d_doppler_th.png',bbox_inches='tight')
plt.show()
sys.exit()








