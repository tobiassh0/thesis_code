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
(nw,nk) = FT2d.shape
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
kernel = 'scharr'
_,kGangle = Kernel(FT2d,kernel=kernel) # scharr or sobel
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
print(kernel+' kernel mean :: ',np.mean(dw_dk))
print(kernel+' kernel medi :: ',np.median(dw_dk))
# plot hist
fig,ax = plt.subplots(figsize=(6,4))
counts,bins,_=ax.hist(dw_dk,bins=1000,range=(-1,1),density=True) # np.log10
dsv = bins[np.argmax(counts)] # doppler shift velocity in units of vA
print(kernel+' kernel max :: ', dsv)
ax.set_ylabel('Normalised count',**tnrfont)
ax.set_xlabel(r'$d\omega/dk$'+'  '+r'$[v_A]$',**tnrfont)
fig.savefig('dw_dk_'+kernel+'_grad.png',bbox_inches='tight')
#plt.show()

plt.imshow((FT2d),**kwargs,extent=[0,kmax/knorm,0,wmax/wnorm])
kx = np.linspace(0,20,100)*knorm
# doppler shifted line 
for i in range(0,int(wmax/wnorm),1):
	w = wnorm*np.ones(len(kx))*i
	ww = w + (dsv*vA)*kx
	plt.plot(kx/knorm,ww/wnorm,color='white',linestyle='--')
plt.xlim(0,20)
plt.ylim(0,10)
plt.ylabel(r'$\omega/\Omega_p$',**tnrfont)
plt.xlabel(r'$kv_A/\Omega_p$',**tnrfont)
plt.plot([0,10],[0,10],color='white',linestyle='--') # vA line
plt.savefig('FT_2d_doppler.png')
#plt.show()
sys.exit()








