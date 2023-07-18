
from func_load import *


def plotSpectroPower(file0,fieldmatrix,times,nx,nt,fs,filefreq=1000,nfft=None,noverlap=None,clim=(None,None),cbar=False):
### In:
#	file0 		: First file in the index list
#	fieldmatrix : Fieldmatrix of a field quantity (typically Bz)
#	times 		: Array of times throughout the simulation
#	nx, nt 		: Number of grid cells and time steps (#files)
#	nfft 		: Length of the FFT [number of files] (default to len(times)//10)
#	noverlap 	: Overlap between FFTs (default this to 6*nfft//8)
#	fs 			: Sampling frequency
#	filefreq 	: The number of files to take spectrograms of (defaults to every 1000th file)
### Out
#	Plots the summed Spectrogram intensity into a power spectra across all frequencies of a signal through time and space

	tc_alpha = 2*const.PI/getCyclotronFreq(file0,'Alpha')
	fc_alpha = getCyclotronFreq(file0,'Alphas')/(2*const.PI)

	if nfft == None:
		nfft = len(times)//10
		if noverlap == None:
			noverlap = 6*nfft//8

	powerheat = np.zeros((nt,nx))
	for i in range(0,nx,filefreq): # loop through
#		if (100*i//nx)%5 == 0: print(100*i//nx,'%') 
		fmx = fieldmatrix[:,i]
		freqs, time, Sxx = signal.spectrogram(fmx,fs,nperseg=nfft,noverlap=noverlap)
		power = []
		for j in range(len(time)):
			power.append(np.sum(Sxx[:,j]**2))#Sxx[freqs, time]
			powerheat[j,i] = np.sum(Sxx[:,j]**2)
		
	L = getGridlen(file0)
	plt.figure(figsize=(8,6))
	powerheat = powerheat[np.all(powerheat!=0,axis=1)]
	im = plt.imshow(powerheat,extent=[0,1,0,times[-1]/tc_alpha],origin='lower',interpolation='nearest',aspect='auto',clim=clim)
	if cbar:
		cbar = plt.colorbar(im)
		cbar.set_label('Spectrogram Power',fontsize=18)
	plt.xlabel(r'$x/L$',fontsize=18)
	plt.ylabel(r'$t/\tau_{c\alpha}$',fontsize=18)

#	cx = [0]
#	x_arr = np.linspace(0,1,100)
#	for c in cx:
#		t_arr = -7*(x_arr-c)+times[-1]/tc_alpha
#		plt.plot(x_arr,t_arr,color='white',linestyle='-.')
#	plt.ylim(0,times[-1]/tc_alpha)
#	plt.xlim(0,1)
	plt.savefig('powerheat.jpeg')
	plt.clf()
	
	return None

def plotStackedSpectro(file0,fieldmatrix,times,nx,nt,fs,freqlim=70,filefreq=1000,pow_t=True,nfft=None,noverlap=None,clim=(None,None),cbar=True,ylims=(0,26E-6)):
### In:
#	file0 		: First file in the index list
#	fieldmatrix : Fieldmatrix of a field quantity (typically Bz)
#	times 		: Array of times throughout the simulation
#	nx, nt 		: Number of grid cells and time steps (#files)
#	nfft 		: Length of the FFT [number of files] (default to len(times)//10)
#	noverlap 	: Overlap between FFTs (default this to nfft//8)
#	fs 			: Sampling frequency
#   freqlim		: Plotting frequency limit (normalised to whatever species you choose)
#	filefreq 	: The number of files to take spectrograms of (defaults to every 1000th file)
### Out
#	Plots the fieldmatrix, Spectrogram intensity and Spectrogram power spectra in a vertically stacked plot. Dependent on the 			position in space

	## choose which power spectra you want to plot
	pow_f = not pow_t

#	tc_alpha = 2*const.PI/getCyclotronFreq(file0,'Alpha')
	tc_alpha = 2*const.PI/getCyclotronFreq(file0,'Deuterons')
	
	fc_D = getCyclotronFreq(sdfread(0),'Deuterons')/(2*const.PI) # hardcoded for now

	if nfft == None: # check to see if parameters are given
		nfft = len(times)//10
		if noverlap == None:
			noverlap = nfft//2

	# initialise Sxmat
	fm_0 = fieldmatrix[:,0]
	freqs, time, Sxx = signal.spectrogram(fm_0,fs,nperseg=nfft,noverlap=noverlap)
	print(freqs,time)
	Sxmat = np.zeros((len(freqs),len(time)))
	for i in range(0,nx,filefreq):
		fmx = fieldmatrix[:,i]
#		plt.figure(figsize=(10,15))
#		plt.subplot(311)
#		plt.plot(times/tc_alpha,(fmx-mean_Bz)/mean_Bz)
#		plt.ylabel(r'$\delta B_z/B_0$',fontsize=18) # change these labels depending on the 
#		plt.xlabel(r'$t/\tau_{c\alpha}$',fontsize=18)
#		ymax = max(np.abs((fmx-mean_Bz)/mean_Bz))
#		plt.ylim(-ymax,ymax) # hardcoded
#		plt.annotate('n={}'.format(i),xy=((times[-1]/tc_alpha)-1.,0.008),xycoords='data',color='k') 
#		##
#		plt.subplot(312)
		freqs, time, Sxx = signal.spectrogram(fmx,fs,nperseg=nfft,noverlap=noverlap)
		Sxmat += Sxx
		print(len(freqs),len(time),Sxx.shape)
		freqs = freqs/fc_D 
#		plt.pcolormesh(time/tc_alpha,freqs,np.log10(Sxx),vmin=clim[0],vmax=clim[1],cmap='Accent')# shading='gouraud'
#		plt.xlim(time[0]/tc_alpha,time[-1]/tc_alpha)
#		plt.ylim(0,70)
#		plt.ylabel(r'$f/f_{cD}$',fontsize=18)
#		plt.xlabel(r'$t/\tau_{c\alpha}$',fontsize=18)
#		if cbar:
#			plt.colorbar()
#		##
#		plt.subplot(313)
#		maxfreqs = []
#		power = []
#		if pow_t:
#			for j in range(len(time)):
##				power.append(np.sum(Sxx[:,j]**2))
#				ind = np.where(Sxx[:,j]**2==max(Sxx[:,j]**2))
#				if freqs[ind[0]]>70: maxfreq = 0
#				else: maxfreq = freqs[ind[0]] 
#				maxfreqs.append(maxfreq)
##			plt.plot(time/tc_alpha,power)
##			plt.ylabel('Spectrogram Power',fontsize=18)
#			plt.plot(time/tc_alpha,maxfreqs)
#			plt.ylabel(r'$f_{max}/f_{cD}$',fontsize=18)
#			plt.xlabel(r'$t/\tau_{c\alpha}$',fontsize=18)
##			plt.ylim(ylims[0],ylims[1])
#			plt.ylim(0,70)
##			plt.savefig('signal_t_{}.jpeg'.format(i))
#			plt.savefig('signal_maxfreq_{}.jpeg'.format(i))
#			plt.clf()
#		if pow_f:
#			for j in range(len(freqs)):
#				power.append(np.sum(Sxx[j,:]**2))
#			plt.plot(power,freqs)
#			plt.xlabel('Spectrogram Power',fontsize=18)
#			plt.ylabel(r'$f/f_{cD}$',fontsize=18)
#			plt.ylim(0,70)
#			plt.xlim(0,1E-23)
#			plt.savefig('signal_f_{}.jpeg'.format(i))
#			plt.clf()
	print('freqs\n',freqs/fc_D)
	dumpfiles([Sxmat,freqs,times],'SxmatCollection_Sft_noverlap_{}_nfft_{}'.format(noverlap,nfft))
#	plt.imshow(np.log10(Sxmat),extent=[0,time[-1]/tc_alpha,0,70],cmap='Accent',interpolation='nearest',origin='lower',aspect='auto')
	plt.pcolormesh(time/tc_alpha,freqs,np.log10(Sxmat),shading='nearest',cmap='jet',clim=(-16.,-9.6))# gouraud
	plt.xlabel(r'$t/\tau_{c\alpha}$',fontsize=18)
	plt.ylabel(r'$f/f_{cD}$',fontsize=18)	
	plt.colorbar()
	plt.xlim(time[0]/tc_alpha,time[-1]/tc_alpha)
	plt.ylim(0,70)
	plt.savefig('sum_Sxmat.png')
	return None

def plotMultipleSx(sim_lst,ratio=True)
	SxMAT = []
	sim_lst = ['traceT_0_00','traceT_0_01','traceT_0_11']
	for sim in sim_lst:
		print(sim)
		sim_loc = getSimulation('/storage/space2/phrmsf/'+sim)
		quant = 'Magnetic_Field_Bz'
		times = read_pkl('times')
		nt = len(times)
		fs = int(len(times)/(times[-1]-times[0]))
		tc_D = 2*const.PI/getCyclotronFreq(sdfread(0),'Deuterons')

		## these are good values
		nfft = len(times)//10 # lower value --> better time resolution (worse freq res)
		noverlap = nfft//2 # better windowing between frames, removes bleeding
		Sxmat, freqs, time = read_pkl('SxmatCollection_Sft_noverlap_{}_nfft_{}'.format(noverlap,nfft))
		freqmax = (70/freqs[-1])*freqs.shape[0]
		SxMAT.append(Sxmat[:int(freqmax),:]) # all time but smaller freq range

	fig,axs = plt.subplots(ncols=3,figsize=(10,6),sharey=True)
	SxMAT=np.array(SxMAT)
	interp = 'nearest'
	cmap = 'jet'

	ratio = False
	if ratio:
		im1 = axs[0].imshow((SxMAT[1,:]/SxMAT[0,:]),origin='lower',interpolation=interp,extent=[0,times[-1]/tc_D,0,70],cmap=cmap,aspect='auto',vmax=3.5,vmin=0.45) # 1%/0%
		im2 = axs[1].imshow((SxMAT[2,:]/SxMAT[0,:]),origin='lower',interpolation=interp,extent=[0,times[-1]/tc_D,0,70],cmap=cmap,aspect='auto',vmax=3.5,vmin=0.45) # 11%/0%
		im3 = axs[2].imshow((SxMAT[2,:]/SxMAT[1,:]),origin='lower',interpolation=interp,extent=[0,times[-1]/tc_D,0,70],cmap=cmap,aspect='auto',vmax=3.5,vmin=0.45) # 11%/1%
		rLabels = [r'$1\%/0\%$',r'$11\%/0\%$',r'$11\%/1\%$']
		cbarLabel = r'$S_{m\%}/S_{n\%}$'
		name = 'SxMat_ratio.png'
	if not ratio:
		im1 = axs[0].imshow(np.log10(SxMAT[0,:]),origin='lower',interpolation=interp,extent=[0,times[-1]/tc_D,0,70],cmap=cmap,aspect='auto',vmax=-9.6,vmin=-16) # 1%/0%
		im2 = axs[1].imshow(np.log10(SxMAT[1,:]),origin='lower',interpolation=interp,extent=[0,times[-1]/tc_D,0,70],cmap=cmap,aspect='auto',vmax=-9.6,vmin=-16) # 11%/0%
		im3 = axs[2].imshow(np.log10(SxMAT[2,:]),origin='lower',interpolation=interp,extent=[0,times[-1]/tc_D,0,70],cmap=cmap,aspect='auto',vmax=-9.6,vmin=-16) # 11%/1%
		rLabels = [r'$0\%$',r'$1\%$',r'$11\%$']
		cbarLabel = r'$\log_{10}(S)$'
		name = 'SxMat_all.png'

	cbar = fig.colorbar(im3)
	cbar.set_label(label=cbarLabel,fontsize=20) 
	for i in range(len(axs)): 
		axs[i].set_title(rLabels[i],fontsize=22)
		axs[i].set_xlabel(r'$t/\tau_{cD}$',fontsize=22)
	axs[0].set_ylabel(r'$f/f_{cD}$',fontsize=22)
	fig.savefig('/storage/space2/phrmsf/'+name)#,bbox_inches='tight')

	return None

def Spectrogram(sim,stacked=True,power=False):
	#### Spectrogram for each k and summation of spectro power
	print(sim)
	sim_loc = getSimulation('/storage/space2/phrmsf/'+sim)
	ind_list = list_sdf(sim_loc)
	quant = 'Magnetic_Field_Bz'
	times = read_pkl('times')
	nt = len(times)
	print(times)
	fs = int(len(times)/(times[-1]-times[0]))
	fc_D = getCyclotronFreq(sdfread(0),'Deuterons')
	tc_D = 2*const.PI/fc_D
	fc_alpha = getCyclotronFreq(sdfread(0),'Alphas')
	tc_alpha = 2*const.PI/fc_alpha #2*const.PI/getCyclotronFreq(sdfread(0),'Alpha')
	## these are good values
	nfft = len(times)//10 # lower value --> better time resolution (worse freq res)
	noverlap = nfft//2 # better windowing between frames, removes bleeding
	maxfreq = 70
	if read:
		Sxmat, freqs, time = read_pkl('SxmatCollection_Sft_noverlap_{}_nfft_{}'.format(noverlap,nfft))
		freqmax = (maxfreq/freqs[-1])*freqs.shape[0]
		Sxmat = Sxmat[:int(freqmax),:] # all time but smaller freq range
		plt.imshow(np.log10(Sxmat),origin='lower',interpolation='nearest',extent=[0,times[-1]/tc_D,0,maxfreq],cmap='jet',aspect='auto',clim=(-9.6,-16))
		plt.colorbar()
		plt.xlabel(r'$t/\tau_{c\alpha}$',fontsize=18)
		plt.ylabel(r'$f/f_{cD}$',fontsize=18)
		plt.savefig('sum_Sxmat.png')	
	else:
		fm = read_pkl('fieldmatrix_'+quant)
		(nt, nx) = fm.shape
		mean_Bz = np.mean(fm[0,:])
		print('fm shape; nx {} , nt {}'.format(nx,nt))
		print('total T :: ',times[-1])
		print('fs :: ', fs)
		print('fc_alpha :: ', fc_alpha)
		print('mean BZ :: ',mean_Bz)
		if stacked:
			plotStackedSpectro(sdfread(0),fm,times,nx,nt,fs,nfft=nfft,noverlap=noverlap,clim=(-22,-11),cbar=True)#,ylims=(0,3E-24))#,clim=(0,5E-13))
		if power:
			plotSpectroPower(sdfread(0),fm,times,nx,nt,fs,filefreq=1,nfft=nfft,noverlap=noverlap,cbar=True)




