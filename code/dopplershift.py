from list_new import *

### General process of doppler shift is a linear trend if the frequency (y-axis) due to the wavenumber (x-axis)
## un-commenting will reveal this graphically with a simple example  
#x = np.arange(0,100,0.01)
#for i in range(0,30,2):
#	y = np.ones(len(x))*i
#	yds = y - x
#	plt.plot(x,yds,color='k',alpha=1-i/30)
#plt.show()

#fig,ax = plt.subplots(figsize=(10,8),ncols=2)
sim_loc = getSimulation('/storage/space2/phrmsf/D_He3_0_10_min_p_0_9')
FT2d = read_pkl('FT_2d_Magnetic_Field_Bz')
d0 = sdfread(0)
klim = 0.5*2*const.PI/getdxyz(d0)
wlim = 0.5*2*const.PI/getdt(d0)
wcp = getCyclotronFreq(d0,'Protons')
LDe = getDebyeLength(d0,'Electrons')
print(wcp,LDe)
#ax[0].imshow(np.log10(FT2d),**kwargs,extent=[0,klim*LDe,0,wlim/wcp],cmap='magma',clim=(-2,6))
#ax[0].set_ylabel(r'$\omega_/\Omega_p$',fontsize=18)
#ax[0].set_yticks([0,5,10,15,20,25,30,35])
#ax[0].set_xlabel(r'$k\lambda_{De}$',fontsize=18)
#ax[0].set_ylim(0,35) ; ax[0].set_xlim(0,0.06)

E0 = 14.57*1e6*const.qe # MeV to Joules
vA = getAlfvenVel(d0)
vp = np.sqrt(2*E0/getMass('Protons'))
vperp = 0.9*vA
vpara = (vp**2 - vperp**2)**.5
print(vA,vp,vA-vp)
ds = vA*((vA-vp)/const.c)
dds = vA*((vA-vpara)/const.c)
kx = np.linspace(0,0.06,100)/LDe
wmax = 20
#for i in range(0,wmax,1):
#	w = wcp*np.ones(len(kx))*i
#	ww = w + ds*kx
#	plt.plot(kx*LDe,ww/wcp,'k--')
##	ww = w + dds*kx
##	plt.plot(kx*LDe,ww/wcp,'b--')
#plt.ylabel(r'$\omega/\Omega_D$',fontsize=18)
#plt.xlabel(r'$k\lambda_{De}$',fontsize=18)
#plt.xlim(0,0.06)
#plt.ylim(0,wmax+5)
#plt.show()

#fig = plt.figure()
#ax = plt.axes(projection='3d')
vds = (vp-vA)*vA/const.c
(NW, NK) = FT2d.shape
print(NW,NK)
v = 0
dv = vds/100
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
#ax[1].imshow(np.log10(tFT2d),**kwargs,extent=[0,klim*LDe,0,wlim/wcp],cmap='magma',clim=(-2,6))
#ax[1].set_yticklabels([])
#ax[1].set_xlabel(r'$k\lambda_{De}$',fontsize=18)
#ax[0].set_ylim(0,35) ; ax[1].set_ylim(0,35)
#ax[0].set_xlim(0,0.06) ; ax[1].set_xlim(0,0.06)
#fig.savefig('/storage/space2/phrmsf/dump/dopplershift_FT2d.png')
powerArr = []
vmax = 1.5*vds
i=0
while v < vmax:
	try:
		powerDop = read_pkl('powerDop_{}'.format(np.around(v,2)))
		powerArr.append(powerDop)
		omegas = np.linspace(0,wmax,len(powerDop))
		varr = np.ones(len(powerDop))*v/vA
#		if i%5==0:
#			ax.plot(omegas,varr,np.log10(powerDop),color='k')
#		else:
#			None
	except:
		print((100*v/vds)//1,'  %')
		tFT2d = np.zeros(FT2d.shape)
		for i in range(tFT2d.shape[1]):
			ki = (i/NK)*klim # real valued
			for j in range(FT2d.shape[0]):
				wj = (j/NW)*wlim # real valued
				wds = wj - ki*v
				jj = (wds/wlim)*NW ## undo doppler shift
				tFT2d[int(jj),i] = FT2d[j,i] # shifted data
		logpowerDop, omegas = powerspectrum(tFT2d,wcp,wlim/wcp,klim*LDe,0,wmax,0,0.06)
		powerDop = 10**logpowerDop
		del logpowerDop
		varr = np.ones(len(powerDop))*v/vA
		dumpfiles(powerDop,'powerDop_{}'.format(np.around(v,2)))
#		ax.plot(omegas/wcp,varr,np.log10(powerDop),color='k')
	v+=dv
	i+=1

powerArr = np.array(powerArr)
im = plt.imshow(np.log10(powerArr),**kwargs,cmap='jet',extent=[0,wmax,0,vmax/vA])
plt.xlabel(r'$\omega/\Omega_p$',fontsize=20)
plt.ylabel(r'$v/v_A$',fontsize=20)
for i in range(0,wmax+1):
	plt.axvline(i,linestyle='--',color='darkgrey')
plt.colorbar(im)
plt.show()
#ax.set_xlabel(r'$\omega/\Omega_p$',fontsize=20)
#ax.set_ylabel(r'$v/v_A$',fontsize=20)
#ax.set_zlabel(r'$\log(p)$',fontsize=20)
#fig.savefig('doppler_power_varr.png',bbox_inches='tight')










