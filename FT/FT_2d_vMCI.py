from func_load import *

sim_lst = ['traceT_0_00','traceT_D_99_T_01','traceT_D_89_T_11','traceT_0_50']

xi2 = [0,1,11,50]
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
	MCI = np.unravel_index(np.nanargmax(FT2d,axis=None), FT2d.shape)
	print(MCI,FT2d[MCI])
	karr=np.linspace(0,kmax,nk) ; warr=np.linspace(0,wmax,nw)
	k_MCI=karr[MCI[1]] # unitless
	w_MCI=warr[MCI[0]] # unitless
#	plt.imshow(np.log10(FT2d),**kwargs,extent=[0,kmax,0,wmax])
#	plt.scatter(k_MCI,w_MCI,marker='s')
#	plt.show()
	# correspond to cold plasma disp
	species = getIonSpecies(sdfread(0))
	omegas = wcyc * np.linspace(0,wmax,100000)
	_,k2,_=coldplasmadispersion(sdfread(0),species[0],species[1],omegas=omegas)
	k2 = k2 * vA/wcyc
	omegas = omegas/wcyc
	kpos = (np.abs(k2 - k_MCI)).argmin()	
	print(k2[kpos])
	plt.scatter(k2-k_MCI,omegas) ; plt.axvline(0,color='darkgrey',linestyle='--') ; plt.show()
	# calculate v_phase (w_MCI / k_MCI)
	v_ph[i] = w_MCI/k_MCI # units of vA
	# calculate gradient (v_group = d w_MCI / d k_MCI)
	grad_inv = np.gradient(k2,omegas[-1]/len(omegas))
	v_gr[i] = 1/grad_inv[kpos] # units of vA
	i+=1

print(v_ph,v_gr)
plt.plot(xi2,v_ph,xi2,v_gr) ; plt.show()
