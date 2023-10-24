from func_load import *

### General process of doppler shift is a linear trend if the frequency (y-axis) due to the wavenumber (x-axis)
## un-commenting will reveal this graphically with a simple example  
#x = np.arange(0,100,0.01)
#for i in range(0,30,2):
#	y = np.ones(len(x))*i
#	yds = y - x
#	plt.plot(x,yds,color='k',alpha=1-i/30)
#plt.show()

def getDopplerShiftICE(sims,FT2darr,normspecies,logthresh=1.8,wkmax=[10,25],plot=True,kernel='custom',dwdkrange=(-2,2),theta=86.3):
	home = os.getcwd()
	# load sim & FT2d 
	fig1,ax1=plt.subplots(nrows=len(sims),figsize=(6,4*len(sims)))
	fig2,ax2=plt.subplots(nrows=len(sims),figsize=(6,4*len(sims)))
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
		tFT2d[FT2d < thresh] = 0
		FT2d = tFT2d # destroy temp arr
		del tFT2d
		
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
		dwdk = np.nan_to_num(1/np.tan(kGangle),posinf=np.nan,neginf=np.nan)
		# remove nan & zero values
		dwdk = dwdk[~np.isnan(dwdk)]
		dwdk = dwdk[dwdk!=0]

		# normalise & thresh
		dw_dk = dwdk * (dw/dk)/vA * (1/10) # TODO: correction factor?
		thresh = (np.abs(dw_dk) < dwdkrange[1])
		dw_dk = dw_dk[thresh]
		print(kernel+' kernel mean :: ',np.mean(dw_dk))
		print(kernel+' kernel medi :: ',np.median(dw_dk))

		# calc histogram
		if plot:
			# plot hist
			counts,bins,_=ax1[i].hist(dw_dk,bins=1000,range=dwdkrange,density=True) # np.log10
			#ax1[i].hist(counts,bins=bins)
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
	fig1.savefig('dw_dk_'+kernel+'_grad.png',bbox_inches='tight')
	ax2[-1].set_xlabel(r'$kv_A/\Omega_p$',**tnrfont)
	fig2.savefig('FT_2d_doppler.png',bbox_inches='tight')
	return None

os.chdir('/storage/space2/phrmsf/lowres_D_He3/')
home = os.getcwd()
sims = [i for i in os.listdir() if 'p_90' in i]
quantity='Magnetic_Field_Bz'
FT2d = []
for sim in sims:
	_=getSimulation(sim)
	FT2d.append(read_pkl('FT_2d_'+quantity))
	os.chdir(home)
getDopplerShiftICE(sims,FT2d,normspecies='Protons')
sys.exit()
ds,va = np.split(dsva,2,axis=1)

