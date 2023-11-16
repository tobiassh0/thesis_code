

def t1():
	wca = 2.1 * (2 * abs(const.qe))/getMass('Alphas')
	w = np.arange(0,40,0.01) * wca
	
	simloc = getSimulation('/storage/space2/phrmsf/traceT/traceT_0_00')
	#dphase,Rtdx = sc.getSpatialCorrelation(sdfread(0),plot=False)
	vA = getAlfvenVel(sdfread(0))
	
	## cold plasma disp frequencies
	_,k2,_= coldplasmadispersion(sdfread(0),w,theta=89.)
	k2 *= (vA/wca)
	plt.xlim(0,100)
	plt.scatter(k2,w/wca,facecolor='b',edgecolor='b')
	
	## 1994 oblique (24)
	theta = 89. * const.PI/180.
	k = np.arange(0,1000,0.01) * (wca/vA)
	kpara = k*np.cos(theta)
	w2 = (k**2 + kpara**2 + (k*kpara)**2 * (vA/wca)**2 )**2 - (2*k*kpara)**2
	tw2 = 0.5*(vA**2)*(k**2 + kpara**2 + (k*kpara)**2 * (vA/wca)**2 + np.sqrt(w2))
	tw = np.sqrt(tw2)
	plt.scatter(k*vA/wca,tw/wca,facecolor='r',edgecolor='r')
	plt.show()

def t2():
	sims = ['A_B_neq','A_B_C_neq','A_B_eq','A_B_C_eq']
	fig,axs = plt.subplots(nrows=2,ncols=2,figsize=(15.5,9),sharex=True,sharey='row')
	fig.subplots_adjust(hspace=0.07,wspace=0.07)
	ax = axs.ravel()
	i=0
	for sim in sims:
		simloc = getSimulation('/home/space/phrmsf/Documents/EPOCH/5_devel/epoch1d/'+sim)
		times = read_pkl('times')
		tcp = 2*const.PI/getCyclotronFreq(sdfread(0),'Protons')
		species = getAllSpecies(sdfread(0))
		print(sim,species)
		if len(species)==3: # protons, protons, electrons
			colors = ['k','r','orange']
			labels = [r'$e$',r'$p$',r'$p$']
			E0 = read_pkl(species[0]+'_KE') # elec
			E1 = read_pkl(species[1]+'_KE') # p
			E2 = read_pkl(species[2]+'_KE') # p
			E3 = np.zeros(len(E0))
			ax[i].plot(times/tcp,E0-np.mean(E0[0:10]),color=colors[0],label=labels[0])
			ax[i].plot(times/tcp,E1-np.mean(E1[0:10]),color=colors[1],label=labels[1])
			ax[i].plot(times/tcp,E2-np.mean(E2[0:10]),color=colors[2],label=labels[2])
			ax[i].legend(loc='best',ncol=len(species),borderpad=0.0,labelspacing=0.05,prop={'size':18})
			ax[i].tick_params(axis='both', which='major', labelsize=18)
		elif len(species)==4: # deuterons, tritons, protons, electrons
			colors = ['k','b','g','orange']
			labels = [r'$e$',r'$D_2$',r'$T_3$',r'$p$']
			E0 = read_pkl(species[0]+'_KE') # elec
			E1 = read_pkl(species[1]+'_KE') # D
			E2 = read_pkl(species[2]+'_KE') # T
			E3 = read_pkl(species[3]+'_KE') # p
			ax[i].plot(times/tcp,E0-np.mean(E0[0:10]),color=colors[0],label=labels[0])
			ax[i].plot(times/tcp,E1-np.mean(E1[0:10]),color=colors[1],label=labels[1])
			ax[i].plot(times/tcp,E2-np.mean(E2[0:10]),color=colors[2],label=labels[2])
			ax[i].plot(times/tcp,E3-np.mean(E3[0:10]),color=colors[3],label=labels[3])
			ax[i].legend(loc='best',ncol=len(species),borderpad=0.0,labelspacing=0.05,prop={'size':18})
			ax[i].tick_params(axis='both', which='major', labelsize=18)
		else: # n/a
			raise SystemExit
		i+=1
	ax[0].set_xlim(0,1)
	ax[0].set_ylabel(r'$\Delta u$'+'  ['+r'$J/m^3$'+']',**tnrfont)
	ax[2].set_ylabel(r'$\Delta u$'+'  ['+r'$J/m^3$'+']',**tnrfont)
	ax[2].set_xlabel(r'$t/\tau_{cp}$',**tnrfont)
	ax[3].set_xlabel(r'$t/\tau_{cp}$',**tnrfont)
	fig.savefig('/home/space/phrmsf/Documents/EPOCH/5_devel/epoch1d/NDW_collection.png',bbox_inches='tight')

def t3():
	theta = np.linspace(0,2*const.PI,100)
	theta_pitch = theta
	TH,TH_P = np.meshgrid(theta,theta_pitch)
	COSDIFF = np.cos(TH-TH_P)
	plt.imshow(COSDIFF,**kwargs,cmap='bwr',extent=[0,2*const.PI,0,2*const.PI])
	plt.show()

def t4():
	# check Delta u_D/Delta u_He3 vs xiHe3
	os.chdir('/storage/space2/phrmsf/lowres_D_He3')
	sims = np.sort([i for i in os.listdir() if 'p_90' in i and i!='0_00_p_90'])
	hlabels = np.array([int(i[2:4]) for i in sims])
	maxduarr,xiarr,marr,qarr=gr.majIons_edens_ratio(sims,species=['Deuterons','He3'],time_norm=r'$\tau_{cD}$',plot=False)
	plt.scatter(xiarr[:,1],maxduarr[:,0]/maxduarr[:,1])
	plt.ylabel('uD/uHe3',**tnrfont)
	plt.xlabel('xiHe3',**tnrfont)
	os.chdir('/storage/space2/phrmsf/lowres_D_He3')
	plt.savefig('maxdu_ratio_vs_xiHe3.png')
	plt.show()

def t5(): # compare growth rates 
	sim_loc = getSimulation('/storage/space2/phrmsf/traceT/old/traceT_highres_0_01')
	theta = 89.
	minions = 'Alphas'
	majions = 'Deuterons'
	wcyca = getCyclotronFreq(sdfread(0),minions)
	omegaall = wcyca*np.linspace(0,25,20000)
	vA = getAlfvenVel(sdfread(0))
	_,kall,_ = coldplasmadispersion(sdfread(0),omegaall)
	emin = 3.5e6
	v0 = np.sqrt(2*emin*const.qe/getMass(minions))
	val = 0.001
	u = 1.67*vA #v0*np.sin(pitch_angle)
	vd = np.sqrt(v0**2 - u**2) #v0*np.sin(pitch_angle)
	vr = v0/1000
	norm_w1, norm_gamma1 = growth_rates_analytical_all(theta, minions, majions , v0, u, kall, omegaall, vr/v0) #v0 in m/s
	w2, gamma2 = growth_rate_man(minions, majions, theta, sdfread(0), u, vd, vr, kall, omegaall)
	w1, gamma1 = norm_w1*wcyca, norm_gamma1*wcyca
	plt.plot(w1/wcyca,gamma1/wcyca,color='r')
	plt.plot(w2/wcyca,gamma2/wcyca,color='b')
	plt.show()
	
if __name__ == '__main__':
	from func_load import *
	#import correlation.spatial_crosscor as sc
	#t1()
	#t2()
	#t3()
	#import energy.gyroresonance as gr
	#t4()
	t5()
