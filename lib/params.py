from list_new import *

#simlst = [,'D_He3_0_10_min_p_0_15']
simlst = ['D_He3_0_10_COLD','D_He3_min_He4_0_8','D_He3_min_He4_0_9','D_He3_0_10_min_p_0_15','D_He3_0_10_min_p_0_9']

for sim in simlst:
	print(sim)
	sim_loc = getSimulation('/storage/space2/phrmsf/'+sim)
	d0 = sdfread(0)
	species = getIonSpecies(d0)
	(xi1, xi2, xi3) = getConcentrationRatios(d0)
	n0 = getMeanquantity(d0,'Derived_Number_Density_Electrons')
	print('n0 : {}'.format(n0))
	print('xi1 : {}\nxi2 : {}\nxi3 : {} : '.format(xi1,xi2,xi3))
	L = getGridlen(d0)
	Nparts = np.zeros(len(species))
	times = read_pkl('times')
	Nt = len(times)
	Nx = len(getGrid(d0)[0])
	print('Nt : {}\nNx : {}'.format(Nt,Nx))
	B0 = getMeanField3D(d0,'Magnetic_Field_B')
	print('B0 : {}'.format(B0))
	for i in range(len(species)):
		try:
			rLspec = getLarmorRadius(d0,species[i])
			print('{} r_L/L: {}'.format(species[i],rLspec/L))	
			Emin = getMeanquantity(d0,'Derived_Average_Particle_Energy_'+species[i])
			print(species[i] + ' E[MeV]: {}'.format(Emin/(1e6*const.qe)))
		except:
			None
		try:
			pw = getQuantity1d(d0,'Particles_Weight_'+species[i])
			Nparts[i] = len(pw)*np.mean(pw) # real No. particles = no. simulated * weight per particle
		except:
			continue # should be 0 for species that arent present
		print('Cspec : '+str(species[i])+' : '+str(Nparts[i]))
	vA = getAlfvenVel(d0)
	print('vA/c : {}'.format(vA/const.c))
	theta = getMagneticAngle(d0)[0]
	print('theta [deg] : {}'.format(theta*180/const.PI))
##	vpara, vperp = getPerpParaVel(d0,'Protons')
##	print('vpara/vA : {}\nvperp/vA : {}'.format(vpara/vA,vperp/vA))
