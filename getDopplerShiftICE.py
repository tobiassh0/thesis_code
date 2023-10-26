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

def getKernelDoppler(sims,FT2darr,normspecies,logthresh=1.8,wkmax=[10,25],plot=True,kernel='custom',dwdkrange=(None,None),theta=86.3):
	home = os.getcwd()
	# load sim & FT2d 
	l = len(sims)
	if l == 1: l+=1
	fig1,ax1=plt.subplots(nrows=l,figsize=(6,4*len(sims)))
	fig2,ax2=plt.subplots(nrows=l,figsize=(6,4*len(sims)))
	#dsva=[] # 2d array of [[normalised gradient, Alfven vel], [... , ...]] per sim
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
		tFT2d = FT2d.copy() # copy arr
		tFT2d[FT2d < logthresh] = 0
		FT2d = tFT2d # replace old FT2d 
		del tFT2d # destroy temp arr
		
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
			counts,bins,_=ax1[i].hist(dw_dk,bins=1000,density=True,range=(-2,2))#dwdkrange) # np.log10
			#ax1[i].hist(counts,bins=bins) dw_dk
			dsv = bins[np.argmax(counts)] # doppler shift velocity, in units of vA
			print(kernel+' kernel max :: ', dsv)
			ax1[i].set_ylabel('Normalised count',**tnrfont)
			#dsva.append([dsv,vA])

			# plot FT2d
			ax2[i].imshow((FT2d),**kwargs,extent=[0,wkmax[1],0,wkmax[0]])
			kx = np.linspace(0,20,100)*knorm
#			kperp = np.sin(theta*const.PI/180)
#			kpara = np.cos(theta*const.PI/180)
#			uperp = 0.9; upara = 6.076 
#			dsth = -(kperp*uperp + kpara*upara)
#			print(dsv,(dsth+1))
			# doppler shifted line 
			for j in range(0,int(wmax/wnorm),1):
				w = wnorm*np.ones(len(kx))*j
				# empirical
				ww = w + (dsv*vA)*kx
				ax2[i].plot(kx/knorm,ww/wnorm,color='white',linestyle='--')
#				# theory
#				tww = w + ((dsth+1)*vA)*kx
#				ax2[i].plot(kx/knorm,tww/wnorm,color='white',linestyle='-.')
			ax2[i].set_xlim(0,20) # reduce plotting limits
			ax2[i].set_ylim(0,10) # " " 
			ax2[i].set_ylabel(r'$\omega/$'+getOmegaLabel(normspecies),**tnrfont)
			ax2[i].plot([0,10],[0,10],color='white',linestyle='--') # vA line
		i+=1
		os.chdir(home)
	print(os.getcwd())
	ax1[-1].set_xlabel(r'$d\omega/dk$'+'  '+r'$[v_A]$',**tnrfont)
	#fig1.savefig('dw_dk_'+kernel+'_grad.png',bbox_inches='tight')
	ax2[-1].set_xlabel(r'$kv_A/\Omega_p$',**tnrfont)
	plt.show()
#	fig2.savefig('FT_2d_doppler.png',bbox_inches='tight')
	return None

#os.chdir('/storage/space2/phrmsf/lowres_D_He3/')
#home = os.getcwd()
#sims = [i for i in os.listdir() if 'p_90' in i]
quantity='Magnetic_Field_Bz'
#FT2d = []
#for sim in sims:
#	_=getSimulation(sim)
#	FT2d.append(read_pkl('FT_2d_'+quantity))
#	os.chdir(home)
#sim=getSimulation('/storage/space2/phrmsf/lowres_D_He3/0_38_p_90')
#sims=[sim]
#FT2d = [read_pkl('FT_2d_'+quantity)]
#getKernelDoppler(sims,FT2d,normspecies='Protons')
#sys.exit()
#ds,va = np.split(dsva,2,axis=1)


def getIntegrateDoppler(sims,FT2darr,normspecies='Protons',wkmax=[20,20],logthresh=1.8):
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
		tsFT2d[sFT2d < logthresh] = 0
		FT2d = tFT2d # replace old FT2d 
		sFT2d = tsFT2d # replace old FT2d 
		del tFT2d, tsFT2d # destroy temp arr

		# extract maximum point along FAW
		argind = np.unravel_index(sFT2d.argmax(),sFT2d.shape)
		ym,xm = argind[0]*swmax/sFT2d.shape[0], argind[1]*skmax/sFT2d.shape[1]
		del sFT2d # remove smaller matrix
		print(xm/knorm,ym/wnorm)
		angles = np.linspace(-const.PI/8,0,10000) # between -45deg and 0deg
		# varying angle of line
		integ = np.zeros(len(angles))
		for i in range(len(angles)):
			# find image edge intercepts
			xi = (xm/knorm-(ym/wnorm)/np.tan(angles[i])) # algebra
			yi = (ym/wnorm-(xm/knorm)*np.tan(angles[i])) # normalised
			# convert to pixel coordinates
			xi*=nk/(wkmax[1]); yi*=nw/(wkmax[0])
			x,y=np.linspace(0,xi,1000),np.linspace(0,yi,1000)
			# map coordinates to integrate
			zi = scipy.ndimage.map_coordinates(FT2d, np.vstack((y,x)))
			# integrate and append to array
			integ[i]=np.sum(zi)
		# maxangle (shared area)
		maxangle = angles[np.argmax(integ)]
		# plot angles vs integral
		fig,ax=plt.subplots(nrows=2,figsize=(6,8))
		ax[0].plot(angles,integ)
		ax[0].plot([0,0],[0,maxangle],color='k',linestyle='--')
		ax[1].imshow(FT2d,**kwargs,extent=[0,wkmax[1],0,wkmax[0]])
		# doppler shifted lines
		dsv = np.tan(maxangle) * (dw/dk)/vA		
		kx = np.linspace(0,wkmax[1],100)*knorm
		for l in range(0,int(wkmax[0]),1):
			w = wnorm*np.ones(len(kx))*l
			ww = w + (dsv*vA)*kx
			ax[1].plot(kx/knorm,ww/wnorm,color='white',linestyle='--')
		
		plt.show()
		sys.exit()
		
			# integrate d,k space along line
		# find maximum power along line as function of angle
		# plot 2d FFT with doppler shifted harmonics
		
sim=getSimulation('/storage/space2/phrmsf/lowres_D_He3/0_38_p_90')
sims=[sim]
FT2d = [read_pkl('FT_2d_'+quantity)]
getIntegrateDoppler(sims,FT2d)
		
		
		