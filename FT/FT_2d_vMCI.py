from func_load import *

sim_lst = ['traceT_0_00','traceT_D_99_T_01','traceT_D_95_T_05','traceT_D_89_T_11','traceT_D_82_T_18','traceT_D_70_T_30','traceT_0_50']
xi2 = np.array([0,0.01,0.05,0.11,0.18,0.30,0.50])

v_ph = np.zeros(len(sim_lst))
v_gr = np.zeros(len(sim_lst))
kmax = 100
wmax = 50
i=0
for sim in sim_lst:
	simloc = getSimulation('/storage/space2/phrmsf/'+sim)
	times = read_pkl('times')
	dt = times[-1]/len(times)
	dx = getdxyz(sdfread(0))
	# load FT2d
	print(sim)
	FT2d = read_pkl('FT_2d_Magnetic_Field_Bz') # assumes already created
	# chopping FT
	klim = 2*0.5*const.PI/dx
	wlim = 2*0.5*const.PI/dt
	vA = getAlfvenVel(sdfread(0))
	wcyc = getCyclotronFreq(sdfread(0),'Deuterons')
	wnorm = wcyc
	knorm = wcyc/vA
	klim_prime = klim/knorm
	wlim_prime = wlim/wnorm
	(nw,nk) = FT2d.shape
	FT2d = FT2d[:int(nw*wmax/wlim_prime),:int(nk*kmax/klim_prime)]
	(nw,nk) = FT2d.shape
	# extract most powerful k and w (MCI wavepacket)
	oset = 5
	MCI = np.unravel_index(np.nanargmax(FT2d[oset:,oset:],axis=None), FT2d[oset:,oset:].shape)
	karr=np.linspace(0,kmax,nk) ; warr=np.linspace(0,wmax,nw)
	k_MCI=karr[MCI[1]+oset] # unitless
	w_MCI=warr[MCI[0]+oset] # unitless
#	plt.imshow(np.log10(FT2d),**kwargs,extent=[0,kmax,0,wmax])
#	plt.scatter(k_MCI,w_MCI,marker='s')
	# correspond to cold plasma disp
	species = getIonSpecies(sdfread(0))
	omegas = wcyc * np.linspace(0,wmax,nk)
	_,k2,_=coldplasmadispersion(sdfread(0),species[0],species[1],omegas=omegas)
	k2 = k2 * vA/wcyc
	omegas = omegas/wcyc
	# calculate v_phase (w_MCI / k_MCI)
	v_ph[i] = w_MCI/k_MCI # units of vA
	# calculate gradient (v_group = [dw/dk]|_MCI)
	dw = (omegas[-1])/len(omegas)
	kgrad = np.gradient(k2,dw)
	kpos = int((np.abs(k2[1:] - k_MCI)).argmin())+1
	print(k2[kpos],k_MCI,w_MCI)
	v_gr[i] = 1/kgrad[kpos] # units of vA
	# correcting vA
#	v_gr[i] = v_gr[i]*vA
#	v_ph[i] = v_ph[i]*vA
	i+=1

## vA lines
#m1 = getMass('Deuterons') ; m2 = getMass('Tritons')
#Z1 = getChargeNum('Deuterons') ; Z2 = getChargeNum('Tritons')
#XI2 = np.linspace(0,1,1000)
#B0 = 2.1 ; n0 = 1e19
#vAarr = B0/np.sqrt(const.mu0*n0*(XI2*(m2-m1*Z2/Z1)+m1/Z1))
#for j in np.arange(0.5,1,0.05):
#	ax.plot(XI2,j*vAarr,linestyle='--',color='darkgrey')
print(v_ph,v_gr)
fig,ax=plt.subplots(figsize=(6,4))
ax.plot(xi2,v_ph,'-s',color='g')
ax.plot(xi2,v_gr,'-s',color='darkorange')
ax.legend(labels=[r'$v_{ph}$',r'$v_{gr}$'],loc='lower right',prop={'size': 18})
ax.set_xlabel(r'$\xi_T$',**tnrfont)
ax.set_ylabel(r'$v/v_A$',**tnrfont)
ax.set_xlim(0,0.51)
ax.set_ylim(0,1)
fig.savefig('/home/space/phrmsf/Documents/thesis_code/ph_gr.png',bbox_inches='tight')
plt.show()

