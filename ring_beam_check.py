
def check_rb(sims,labels,minspec='Protons'):
	frac = 1000
	fig,axs = plt.subplots(nrows=3,ncols=len(sims),figsize=(10,10))
	figg = plt.figure()
	ax = figg.add_subplot(projection='3d')
	for j in range(len(sims)):
		simloc = getSimulation(sims[j])
		d0 = sdfread(0)
		p_x = getQuantity1d(d0,'Particles_Px_'+minspec)
		p_y = getQuantity1d(d0,'Particles_Py_'+minspec)
		p_z = getQuantity1d(d0,'Particles_Pz_'+minspec)
		ax.scatter(p_x,p_y,p_z)
		p_arr = np.array([p_x[::frac],p_y[::frac],p_z[::frac]])
		dim = ['x','y','z']
		for i in range(0,3):
			axs[i,j].hist(p_arr[i],bins=500)
			axs[i,j].set_title(str(labels[j])+' '+dim[i])
		os.chdir('..')
		print(labels[j])
		figg.savefig('scatter_p_{}.png'.format(sims[j]))
		figg.clf()
	# fig.savefig('hist_p.png')

	return fig

if __name__=='__main__':
	from func_load import *
	home = os.getcwd()
	# # D-He3
	# os.chdir('/storage/space2/phrmsf/lowres_D_He3')
	# sims = np.sort([i for i in os.listdir() if 'p_90' in i])
	# hlabels = np.array([int(i[2:4]) for i in sims])
	# check_rb(sims,labels=hlabels)
	# # D-T
	# os.chdir('/storage/space2/phrmsf/traceT')
	# sims = [i for i in os.listdir() if 'traceT' in i]
	# sims.sort()
	# sims.append(sims[0]) # add 0% onto end
	# sims = np.array(sims[1:]) # remove duplicate and make numpy array
	# tlabels = np.array([int(i[-2:]) for i in sims])
	# thresh = tlabels >= 18
	# check_rb(sims[thresh],labels=tlabels[thresh],minspec='Alphas')
	# Fast electrons
	os.chdir('/storage/space2/phrmsf/ECRH')
	sims = np.sort([i for i in os.listdir() if 'ECRH_JT60U_' in i and 'images' not in i])
	labels = np.array([i[-1:] for i in sims])
	fig = check_rb(sims,labels,minspec='FElectrons')
	os.chdir(home)
	fig.savefig('hist_p.png')

