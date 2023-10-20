from func_load import *

home = '/storage/space2/phrmsf/lowres_D_He3/'
sim_lst = ['0_15_p_90','0_22_p_90','0_25_p_90']

for sim in sim_lst:
	simloc = getSimulation(home+sim)
	ind_lst = list_sdf(simloc)
	times = read_pkl('times')
	troll = np.roll(times,1)
	dtarr = times-troll
	dtarr = dtarr[1:]
	plt.plot(ind_lst[1:],dtarr)
	plt.ylim(1.46e-11,1.5e-11)
	os.chdir('/home/space/phrmsf/Documents/thesis_code/')
	plt.savefig(sim+'_dtarr.png')
	print(sim+'  fig saved.')
	plt.clf()
	
