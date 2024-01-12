
def funcmean(data):
	return np.mean(data)

def funcvar(data):
	return np.var(data)

def funcskew(data):
	return stats.skew(data)

def funckurtosis(data):
	return stats.kurtosis(data)

def getMoment(name,data):
	mom_dict = {
	"mean" : funcmean,
	"var" : funcvar,
	"skw" : funcskew,
	"kur" : funckurtosis
	}
	moment = mom_dict.get(name)
	return moment(data)

def velocity_bulk(sim,restart_files=[],dim=['x','y','z'],majspec='Deuterons',moments=['mean','var','skw','kur']):
	# nrows = dimensions, ncols = restart files
	colors=['r','g','b','k']
	fig,axs=plt.subplots(ncols=len(dim),figsize=(12,8))
	axs.ravel()
	simloc=getSimulation(sim)
	vth = getThermalVelocity(sdfread(0),majspec)
	pnorm = getMass(majspec)*vth
	for i in range(len(restart_files)):
		di = sdfread(restart_files[i])
		for m in range(len(moments)):
			for j in range(len(dim)):
				pmin = getQuantity1d(di,'Particles_P'+dim[j]+'_'+majspec)/pnorm
				moment = getMoment(moments[m],pmin)
				axs[j].scatter(restart_files[i],moment,color=colors[m])
			print(restart_files[i])
	# plt.xlim(-6,6)
	fig.savefig(sim+'_mom_rest.png')
	plt.show()
	return None


if __name__=='__main__':
	from func_load import *
	vel = 'Particles_P'
	dim = ['x','y','z']
	# D-He3
	os.chdir('/storage/space2/phrmsf/lowres_D_He3/')
	home=os.getcwd()
	sims = ['0_00_p_90']
	majspec='Deuterons'
	# read restart.visit files
	for sim in sims:
		simloc=getSimulation(sim)
		with open('restart.visit','r') as file:
			# read file, replace sdf with '', split str into list, remove empties (list(filter(...))
			trestart_files=list(filter(None,file.read().replace('.sdf','').split('\n')))
			# remove TempOut files if used
			for f in trestart_files:
				if 'TempOut' in f: trestart_files.remove(f)
			restart_files = [int(i) for i in trestart_files] # convert to int
			del trestart_files
			file.close()
		os.chdir(home)
		velocity_bulk(sim,restart_files)

