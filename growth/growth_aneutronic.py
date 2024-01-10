
def growth_aneutronic(sims,labels):
	ls = len(sims)
	fig,axs=plt.subplots(nrows=3,ncols=ls//3,figsize=(8,5))
	# extract basic parameters (how to include second species?)
	for i in range(ls):
		simloc = getSimulation(sims[i])
		d0 = sdfread(0)
		species = getIonSpecies(d0)
		majions, maj2ions, minions = species # split species 
		theta_deg,_ = getMagneticAngle(d0)
		xi1,xi2,xi3 = getConcentrationRatios(d0)
		vA = getAlfvenVel(d0)
		print(xi2/xi1,vA/const.c)

#		print(((100*xi2)//1)/100) # round to nearest 100th
#		posomega, posgamma = growth_rate_manual(minions=minions,majions=majions,maj2ions=maj2ions,wmax=wmax,theta_deg=theta_deg,xi3=10**(-4),xi2=)
#		axs[i].plot(posomega,posgamma)

		os.chdir('..') # change back to home dir 

	return None

if __name__=='__main__':
	from func_load import * 
	os.chdir('/storage/space2/phrmsf/lowres_D_He3/')
	sims = np.sort([i for i in os.listdir() if 'p_90' in i])
	hlabels = [int(i[2:4]) for i in sims]
	print(sims,hlabels)
	growth_aneutronic(sims,hlabels)
