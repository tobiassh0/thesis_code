from func_load import *
import scipy
### General process of doppler shift is a linear trend if the frequency (y-axis) due to the wavenumber (x-axis)
## un-commenting will reveal this graphically with a simple example  
#x = np.arange(0,100,0.01)
#for i in range(0,30,2):
#	y = np.ones(len(x))*i
#	yds = y - x
#	plt.plot(x,yds,color='k',alpha=1-i/30)
#plt.show()

def getKernelDoppler(sims,FT2darr,normspecies,wkmax=[10,25],logthresh=1.8,kernel='custom',dwdkrange=(-0.5,0.5),theta=86.3,plot=True):
	"""
	Uses the kernel method (see list_new.py and this link [https://pyimagesearch.com/2021/05/12/image-gradients-with-opencv-sobel-and-scharr/]) 
	to find the gradients within the image. Uses a custom/modified Scharr kernel which is more sensitive to 
	-45deg gradients. Plots the histogram of these gradients as well as the converted grad to vA atop the 2d 
	FFT.
		params in
			sims:			list of sims to analyse and extract each gradient from
			FT2darr:		list of 2d FFTs corresponding to the sims in list sims 
			normspecies:the normalisation species used for w,k space
			wkmax:		the limits of the 2d FFT [wmax,kmax] which to plot, in units of normspecies normalisation 
			logthresh:	the log value of the thresh which to apply to the 2d FFT 
			kernel:		the kernel to use (sobel, scharr, custom)
			dwdkrange:	the range over dw/dk (grad) space which to calculate the histogram (modify this for accurate working)
			theta:		redundant parameter, but useful if you want to plot the theoretical doppler as well
		params out
			Plots and saves the histograms of gradients (dw/dk) and 2d FFT with doppler shifted harmonics for the whole 
			array of sims provided. 
			dsvarr:		an array of all doppler velocities in units of vA (per sim)
	"""	
	home = os.getcwd()
	# load sim & FT2d 
	l = len(sims)
	if l == 1: l+=1
	fig1,ax1=plt.subplots(nrows=l,figsize=(6,4*len(sims)))
	fig2,ax2=plt.subplots(nrows=l,figsize=(6,4*len(sims)))
	dsvarr=[] # 1d array of gradient in units of m/s per sim
	i=0
	for sim in sims:
		## calc gradients in image
		# setup
		sim_loc = getSimulation(sim)	
		d0 = sdfread(0)
		times = read_pkl('times')
		vA = getAlfvenVel(d0)
		print(vA/const.c)
		klim = 0.5*2*const.PI/getdxyz(d0)
		wlim = 0.5*2*const.PI/getdt(times)
		wnorm = getCyclotronFreq(d0,normspecies)
		knorm = wnorm/vA
		# freq resolutions (no factor half)
		dk = 2*const.PI/getGridlen(d0)
		dw = 2*const.PI/times[-1]
	
		# cut FT2d into size needed
		FT2d = FT2darr[i]
		(nw,nk) = FT2d.shape
		wmax = wkmax[0]*wnorm ; kmax = wkmax[1]*knorm
		FT2d = np.log10(FT2d[:int(nw*wmax/wlim),:int(nk*kmax/klim)])
		(nw,nk) = FT2d.shape
		
		# threshold FT2d
		tFT2d= FT2d.copy() # copy arr
		tFT2d[FT2d < logthresh] = 0
		ttFT2d = tFT2d.copy()
		ttFT2d[FT2d > logthresh] = 1
		FT2d = ttFT2d # replace old FT2d 
		del ttFT2d # destroy temp arr
		
		# Kernel gradient map
		_,kGangle = Kernel(FT2d,kernel=kernel) # scharr or sobel

		# gradients as angles
		kGangle = kGangle[1:-1,1:-1]# remove abberations around edge
#		plt.imshow(kGangle,extent=[0,wkmax[1],0,wkmax[0]],**kwargs)	; plt.show()	
		kGangle = kGangle.flatten()

		# convert to all negative angles (easier to calc real gradient)
		for g in range(len(kGangle)):
			if kGangle[g] > 0:
				kGangle[g]-=const.PI #rad

		# convert grad from angle to velocity (convert +-inf to nan)
		dwdk = np.nan_to_num(np.tan(kGangle),posinf=np.nan,neginf=np.nan)
		# remove nan & zero values
		dwdk = dwdk[~np.isnan(dwdk)]
		dwdk = dwdk[dwdk!=0]

		# normalise & thresh
		dw_dk = dwdk * (dw/dk)/vA # TODO: correction factor?
		#thresh = (np.abs(dw_dk) < dwdkrange[1])
		#dw_dk = dw_dk[thresh]
		print(kernel+' kernel mean :: ',np.mean(dw_dk))
		print(kernel+' kernel medi :: ',np.median(dw_dk))

		# calc histogram
		if plot:
			# plot hist
			counts,bins,_=ax1[i].hist(dw_dk,bins=1000,density=True,range=dwdkrange) # np.log10
			#ax1[i].hist(counts,bins=bins) dw_dk
			dsv = bins[np.argmax(counts)] # doppler shift velocity, in units of vA
			dsvarr.append(dsv*vA)
			print(kernel+' kernel max :: ', dsv)
			ax1[i].set_ylabel('Normalised count',**tnrfont)
			#dsva.append([dsv,vA])

			# plot FT2d
			ax2[i].imshow((tFT2d),**kwargs,extent=[0,wkmax[1],0,wkmax[0]])
			del tFT2d
			kx = np.linspace(0,20,100)*knorm
			kperp = np.sin(theta*const.PI/180)
			uperp = 0.9; upara = 6.076 
			dsth = -(kperp*uperp)# + kpara*upara)
#			print(dsv,(dsth+1))
			# doppler shifted line 
			for j in range(0,int(wmax/wnorm),1):
				w = wnorm*np.ones(len(kx))*j
				# empirical
				ww = w - (abs(dsv)*vA)*kx
				ax2[i].plot(kx/knorm,ww/wnorm,color='white',linestyle='--')
#				# theory
#				tww = w - (abs(dsth+1)*vA)*kx
#				ax2[i].plot(kx/knorm,tww/wnorm,color='white',linestyle='-.')
			ax2[i].set_xlim(0,20) # reduce plotting limits
			ax2[i].set_ylim(0,10) # " " 
			ax2[i].set_ylabel(r'$\omega/$'+getOmegaLabel(normspecies),**tnrfont)
			ax2[i].plot([0,10],[0,10],color='white',linestyle='--') # vA line
		i+=1
		os.chdir(home)
	print(os.getcwd())
	ax1[-1].set_xlabel(r'$d\omega/dk$'+'  '+r'$[v_A]$',**tnrfont)
	fig1.savefig('dw_dk_'+kernel+'_grad.png',bbox_inches='tight')
	ax2[-1].set_xlabel(r'$kv_A/\Omega_p$',**tnrfont)
#	plt.show()
	fig2.savefig('FT_2d_doppler.png',bbox_inches='tight')
	return dsvarr


def getIntegrateDoppler(sims,FT2darr,normspecies='Protons',wkmax=[20,20],logthresh=1.8):
	"""
	Extracts the maximum point in a smaller array (so don't use edge of array and bias sample) of the 2d FFT
	the loops through multiple angles (angle +ve clockwise from North) and extracts the values of the threshed
	2d FFT array, summates them then plots this integrand vs its corresponding angle. Can then find the angle
	of doppler shift and convert this to a velocity, plotting this atop the 2d FFT for multiple l 
	harmonics.  
		params in
			sims:			list of sims to analyse and extract each gradient from
			FT2darr:		list of 2d FFTs corresponding to the sims in list sims 
			normspecies:the normalisation species used for w,k space
			wkmax:		the limits of the 2d FFT [wmax,kmax] which to plot, in units of normspecies normalisation 
			logthresh:	the log value of the thresh which to apply to the 2d FFT 
		params out
			Plots the gradient angle vs. normalised integral to the number of cells np.sum(zi)/len(zi) 
	"""	
	for i in range(len(sims)):
		# setup
		sim_loc = getSimulation(sim)	
		d0 = sdfread(0)
		times = read_pkl('times')
		vA = getAlfvenVel(d0)
		print(vA/const.c)
		klim = 0.5*2*const.PI/getdxyz(d0)
		wlim = 0.5*2*const.PI/getdt(times)
		wnorm = getCyclotronFreq(d0,normspecies)
		knorm = wnorm/vA
		# freq resolutions (no factor half)
		dk = 2*const.PI/getGridlen(d0)
		dw = 2*const.PI/times[-1]
	
		# cut FT2d into size needed
		FT2d = FT2darr[i]
		(nw,nk) = FT2d.shape
		wmax = wkmax[0]*wnorm ; kmax = wkmax[1]*knorm
		swmax = (wkmax[0]-10)*wnorm ; skmax = (wkmax[1]-10)*knorm
		# reduce limit and find max
		sFT2d = np.log10(FT2d[:int(nw*swmax/wlim),:int(nk*skmax/klim)])
		FT2d = np.log10(FT2d[:int(nw*wmax/wlim),:int(nk*kmax/klim)])
		(nw,nk) = FT2d.shape

		# threshold FT2d
		tFT2d = FT2d.copy() ; tsFT2d = sFT2d.copy() # copy arr
		tFT2d[FT2d < logthresh] = 0
		tFT2d[FT2d > logthresh] = 1
		tsFT2d[sFT2d < logthresh] = 0
		tsFT2d[sFT2d > logthresh] = 1
		FT2d = tFT2d # replace old FT2d 
		sFT2d = tsFT2d # replace old FT2d 
		del tFT2d, tsFT2d # destroy temp arr

		# extract maximum point along FAW
		argind = np.unravel_index(sFT2d.argmax(),sFT2d.shape)
		ym,xm = argind[0]*swmax/sFT2d.shape[0], argind[1]*skmax/sFT2d.shape[1]
		del sFT2d # remove smaller matrix
		print(xm/knorm,ym/wnorm)
		angles = np.linspace(const.PI/2,const.PI,1000) # between -180deg and -90deg
		# varying angle of line
		integ = np.zeros(len(angles))
		for i in range(len(angles)):
			# find image edge intercepts
			xi = (xm/knorm-(ym/wnorm)/np.tan(angles[i])) # algebra
			yi = (ym/wnorm-(xm/knorm)*np.tan(angles[i])) # normalised
#			ax.plot([0,xi],[yi,0],linestyle='--',color='white')
			# convert to pixel coordinates
			xi*=nk/(wkmax[1]); yi*=nw/(wkmax[0])
			x,y=np.linspace(0,xi,1000),np.linspace(0,yi,1000)
			# map coordinates to integrate
			zi = scipy.ndimage.map_coordinates(FT2d, np.vstack((y,x)))
			# integrate and append to array
			integ[i]=np.sum(zi/len(zi)) # normalise to integ per cell
		# maxangle (shared area)
		maxangle = angles[np.argmax(integ)]
		# plot angles vs integral
		fig,ax=plt.subplots(nrows=2,figsize=(6,8))
		ax[0].plot(angles*180/const.PI,integ)
		ax[0].plot([0,0],[0,maxangle],color='k',linestyle='--')
		ax[1].imshow(FT2d,**kwargs,extent=[0,wkmax[1],0,wkmax[0]])
		# doppler shifted lines
		dsv = np.tan(maxangle) * (dw/dk)/vA		
		kx = np.linspace(0,wkmax[1],100)*knorm
		for l in range(0,int(wkmax[0]),1):
			w = wnorm*np.ones(len(kx))*l
			ww = w + (dsv*vA)*kx
			ax[1].plot(kx/knorm,ww/wnorm,color='white',linestyle='--')
		ax[1].set_xlabel(r'$kv_A/$'+getOmegaLabel(normspecies),**tnrfont)
		ax[1].set_ylabel(r'$\omega/$'+getOmegaLabel(normspecies),**tnrfont)
		ax[1].set_ylim(0,wkmax[0]) ; ax[1].set_xlim(0,wkmax[1])

		plt.show()
		sys.exit()
	return None	
		
quantity='Magnetic_Field_Bz'

## kernel
#os.chdir('/storage/space2/phrmsf/lowres_D_He3/')
#home = os.getcwd()
#sims = np.sort([i for i in os.listdir() if 'p_90' in i])
#FT2d = []
#for sim in sims:
#	_=getSimulation(sim)
#	FT2d.append(read_pkl('FT_2d_'+quantity))
#	os.chdir(home)
#dsvarr = getKernelDoppler(sims,FT2d,normspecies='Protons')
#sys.exit()

## line integrate
sim=getSimulation('/storage/space2/phrmsf/lowres_D_He3/0_38_p_90')
sims=[sim]
FT2d = [read_pkl('FT_2d_'+quantity)]
getIntegrateDoppler(sims,FT2d)
		
		
		