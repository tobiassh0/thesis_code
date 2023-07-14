import numpy as np
import matplotlib.pyplot as plt 
import my_constants as const
from list_new import *
from scipy import signal
from scipy.fft import fftshift
import time
from scipy.optimize import curve_fit

def mvavg(arr,points=np.linspace(-3,3,7,dtype=int)):
	# IN 
	#	arr : 1d array
	#	points : list of points to move forward and backwards 
	#			across the data to take an average of. 
	# OUT
	#	new_arr : The moving average array
	n = len(arr)
	new_arr = np.zeros(n)
	for i in range(0,n):
		for p in points:
			if i+p > n-1 or i+p < 0:
				new_arr[i] = arr[i]
			else:
				new_arr[i] += arr[i+p]
		new_arr[i]=new_arr[i]/len(points)

	return new_arr

def tophat(length, width, offset=0, height=1):
	# IN
	#	length : the length of the top-hat array (including 0s)
	# 	width : the wdith of the top-hat for when it's >0 
	# 	offset : when to begin the top-hat in your data series (defaults to 0)
	#	height : the non-zero value of the top-hat array (defaults to 1)
	# OUT
	#	toph : the array of length with 0s and height for some width in your signal
	toph = []
	for i in range(length):
		if i/length > offset and i/length < offset+width:
			toph.append(height)
		else:
			toph.append(0)
	return toph

def mask_area(data, xlim, ylim, left, bottom, width, height, mval=0, cbar=False):
	# Mask an area within the data with a square shape starting at bottom left 
	# and moving right and upwards by width and height (in units of the axes)
	# Replaces the data array values within the shape with mval.
	# IN
	#	data : 2d array of data with shape nx, ny
	#	xlim, ylim : limits of the data in x and y (in units of the axes)
	#	left, bottom, width, height : limits of the masking area (box)
	#	mval : value to replace data within the mask
	# OUT
	#	fig : figure of the new data with an applied mask
	
	(nx, ny) = data.shape
	mask_x = nx * (left/xlim)
	mask_y = ny * (bottom/ylim)
	mask_x_width = nx * (width/xlim)
	mask_y_width = ny * (height/ylim)
	fig, axs = plt.subplot(1,2)
	ax = axs[0]
	ax.imshow(data, extent=[0,xlim,0,ylim], interpolation='nearest', origin='lower', cmap='magma')
	ax = axs[1]
	data[int(mask_x):int(mask_x+mask_x_width),int(mask_y):int(mask_y+mask_y_width)] = mval
	ax.imshow(data, extent=[0,xlim,0,ylim], interpolation='nearest', origin='lower', cmap='magma')
	if cbar:
		plt.colorbar()
	return fig

##################################
mD = const.me_to_mD
mT = const.me_to_mT

##################################

JETdata = np.loadtxt('JET26148_ICE_POWER.txt',delimiter=',')
JETpower, JETfreqs = JETdata[:,1], JETdata[:,0] # 2 columns, N rows
wnorm = 2*const.PI*17E6
print('wnorm [MHz] :: {}'.format(wnorm/(2*const.PI*1e6)))
JETfreqs = 2*const.PI*JETfreqs*1e6/(wnorm) # convert MHz to wcD
maxJETfreqs = round(max(JETfreqs))
print('MAX FREQ (JET wcD):: ',maxJETfreqs)
fig,ax=plt.subplots(figsize=(8,4))
ax.annotate(r'$\Omega_D=17$'+'MHz',xy=(0.1,0.8),xycoords='axes fraction',fontsize=20)
ax.plot(JETfreqs,JETpower,'k-')
ax.set_xlabel(r'$\omega/\Omega_D$',fontsize=20)
ax.set_ylabel('ICE intensity [dB]',fontsize=20)
ax.set_xlim(0,maxJETfreqs)
fig.savefig('JET26148_power.png',bbox_inches='tight')
fig.savefig('JET26148_power.eps',bbox_inches='tight')

##################################
##### Plot 2d FFT 
#### load sim
#sim_lst = ['traceT_0_11']#,'traceT_0_01','traceT_0_00','cold_JET26148']
#for sim in sim_lst:
#	sim_loc = getSimulation('/storage/space2/phrmsf/'+sim)
#	FT_2d = read_pkl('FT_2d_Magnetic_Field_Bz')
#	times= read_pkl('times')
#	(nw, nk) = FT_2d.shape
#	print(nw,nk)
#	nx, nt = 2*nk, len(times)
#	d0 = sdfread(0)
#	T = times[-1]
#	L = getGridlen(d0)
#	### get consts
#	species = getIonSpecies(d0)
#	wc_maj = getCyclotronFreq(d0,species[0])
#	wce = getPlasmaFreq(d0,'Electrons')
#	wpi = getPlasmaFreq(d0, species[0])
#	wpe = getPlasmaFreq(d0, 'Electrons')
#	lambdaD = getDebyeLength(d0,'Electrons')
#	va = getAlfvenVel(d0)
#	wnorm = wc_maj
#	knorm = wnorm/va
#	dx = getdxyz(sdfread(0))
#	print('dx : {}'.format(dx))
#	### normalise lims
#	klim, wlim = 0.5*nx*2*const.PI/L, 0.5*nt*2*const.PI/T 
#	klim_prime, wlim_prime = klim/knorm, wlim/wnorm
#	in_klimprime = 100
#	in_wlimprime = 50
#	### cut FT_2d to area of interest
#	w_lim, k_lim = FT_2d.shape[0]*(in_wlimprime/wlim_prime), FT_2d.shape[1]*(in_klimprime/klim_prime)
#	FT_2d = FT_2d[:int(w_lim),:int(k_lim)]
#	fig, ax = plot2dTransform(FT_2d, [va,True], in_klimprime, in_wlimprime, Omega_label=r'$\Omega_D$' ,cbar=False, clim=(-3,6), cmap='magma')
#	#k = np.linspace(0,100,1000)
#	#khatch = [0,38,42,45.5,49]
#	#m = 0.45
#	#for kdp in khatch:
#	#	ax.plot(k,m*(k-kdp),color='k')
#	#ax.plot(k,k,color='k',linestyle='--')
#	#### plot horizontal harmonics
#	#mT = const.me_to_mT
#	#mD = const.me_to_D2
#	#for i in range(0,10):
#	#	plt.axhline(i,alpha=0.5,color='w',linestyle='--')
#	#	plt.axhline(i*mT/mD,alpha=0.5,color='w',linestyle='-.')
#	### Cold plasma dispersion
#	omegas = wnorm*np.linspace(0,in_wlimprime,10000) # do this so that the range of k near w=0 does not diverge, then shift it down
#	k1,k2,k3=coldplasmadispersion(d0,'Deuterons','Electrons',z1=1,z2=1,omegas=omegas) # two solutions to the cold plasma dispersion
#	knorm = (va/wc_maj)
#	k1,k2,k3 = k1*knorm, k2*knorm, k3*knorm
#	print(k2)
#	thresh = k2 > 0 # threshold the array so it only plots the FAW and not the horizontal line to the 0 parts of the dispersion
#	ax.plot(k2[thresh],omegas[thresh]/wnorm,color='k',linestyle='-',alpha=0.75)
#	print('k1 : ',k1,'\n','k2 : ',k2)
#	ax.set_ylim(0,in_wlimprime)
#	ax.set_xlim(0,in_klimprime)
#	### Wavemodes
#	ax = ColdWaveModes(ax,[wpi, wpe, wc_maj, wce],va,wnorm,LHact=True,W2=True) ## check list_new to see which modes we can plot (bool)
#	#plt.show()
#	plotting(fig,ax,'FT_2d_Magnetic_Field_Bz')

###############################################################################################################
###############################################################################################################

#### Waiting for sim to finish
#sim_lst = ['BoxTest']
#for sim in sim_lst:
#	print(sim)
#	sim_loc = getSimulation('/home/space/phrmsf/Documents/EPOCH/5_devel/epoch1d/'+sim)
#	short = True
#	print('Running FT_2d, energy and power :\nWaiting...')
#	TotFileNum = 2001
#	while short: # continuously reading file dir, will execute when file threshold (hardcoded) is reached
#		try:
#			ind_list = list_sdf(sim_loc)
#		except:
##			print('lock test file found')
#			ind_lst = [0]
#			time.sleep(100) # wait 100 seconds before trying again
#		if len(ind_list) < TotFileNum:
#			time.sleep(5*60) # wait 5 minutes before trying again
#			continue			
#		else:	
#			print('Executing...')
#			short = False
#	quant = 'Magnetic_Field_Bz'
#	file0 = sdfread(0)
#	spec_list = getIonSpecies(file0)
#	Nx = len(getGrid(file0)[0])
#	maj_species = spec_list[0]
#	batch_size = 500
#	fieldmatrix, times = load_tot_fieldmatrix(ind_list,Nx,file0,quant,batch_size)
#	dumpfiles(fieldmatrix, 'fieldmatrix_'+quant)
#	dumpfiles(times, 'times')
#	fieldmatrix = read_pkl('fieldmatrix_'+quant)
#	FT_2d = get2dTransform(fieldmatrix,window=True)
#	dumpfiles(FT_2d,'FT_2d_'+quant)
#	times = read_pkl('times')

#	dlast = sdfread(ind_list[-1])
#	wc_maj = getCyclotronFreq(file0,'Deuterons')
#	wce = getPlasmaFreq(file0,'Electrons')
#	wpi = getPlasmaFreq(file0, 'Deuterons')
#	wpe = getPlasmaFreq(file0, 'Electrons')
#	wnorm = wc_maj
#	va = getAlfvenVel(file0)
#	lambdaD = getDebyeLength(file0,'Electrons')
#	klim, wlim = batch_getDispersionlimits((ind_list,file0,dlast,times)) # non-normalised units
#	DISP_DATA = ind_list, file0, dlast, times, klim, wlim, wc_maj, wce, va, lambdaD, wpe, wpi, wnorm
#	klim_prime, wlim_prime, tlim_prime = norm_DispersionLimits(DISP_DATA)#,species=maj_species)
#	
#	in_klimprime = 100
#	in_wlimprime = 50
#	
#	w_lim, k_lim = FT_2d.shape[0]*(in_wlimprime/wlim_prime), FT_2d.shape[1]*(in_klimprime/klim_prime)
#	FT_2d = FT_2d[:int(w_lim),:int(k_lim)]
#	print('plotting shape: ',np.shape(FT_2d))
#	fig, ax = plot2dTransform(FT_2d, va_wci=[va,wnorm==wc_maj], klim=in_klimprime, wlim=in_wlimprime, Omega_label=getOmegaLabel(maj_species),cmap='magma')
#	plotting(fig,ax,'FT_2d_'+quant)	

#	energy_int = energies(sim_loc,min_species=spec_list[2],maj_species=spec_list[0],maj2_species=spec_list[1],frac=1,plot=True,integ=True)
#	power(klim_prime=klim_prime,wlim_prime=wlim_prime,wmax=35,kmax=in_klimprime,norm_omega=getOmegaLabel(maj_species),quant=quant,plot=True)
#	

###############################################################################################################
###############################################################################################################

#sim_lst = ['BoxTest','traceT_D_89_T_11'] # do Box test separate
#fig,axs = plt.subplots()
#fig.subplots_adjust(hspace=0.28)
#ax1 = plt.subplot2grid((2, 2), (0, 0),colspan=2) # setting colspan=2 will
#ax2 = plt.subplot2grid((2, 2), (1, 0))           # move top pie chart to the middle
#ax3 = plt.subplot2grid((2, 2), (1, 1), sharey=ax2)

### small BoxTest
#sim_loc = getSimulation('/home/space/phrmsf/Documents/EPOCH/5_devel/epoch1d/BoxTest')
#fm=read_pkl('fieldmatrix_Magnetic_Field_Bz')
#times = read_pkl('times')
#tc_D = 2*const.PI/getCyclotronFreq(sdfread(0),'Deuterons')
#(nt,nx) = fm.shape
#print(nt,nx)
#Lmin = getGridlen(sdfread(0))
#Tmin = times[-1]/tc_D
#print('Lmin : {}, Tmin : {}'.format(Lmin,Tmin))
#ax3.imshow(fm,extent=[0,Lmin,0,Tmin],origin='lower',interpolation='nearest',aspect='auto',cmap='jet')
#xpos=Lmin/10 ; ypos=Tmin/10
#ax3.set_xlabel(r'$x$'+'  '+'[m]',fontsize=18)
##ax3.set_ylabel(r'$t/\tau_{cD}$',fontsize=18)
#ax3.text(xpos,ypos,'(c)',color='k',fontsize=18)

### large sim, full
#sim_loc = getSimulation('/storage/space2/phrmsf/'+sim_lst[1])
#fm=load_batch_fieldmatrix(quantity='Magnetic_Field_Bz',para=False) #read_pkl('fieldmatrix_Magnetic_Field_Bz')
#times = read_pkl('times')
#tc_D = 2*const.PI/getCyclotronFreq(sdfread(0),'Deuterons')
#(nt,nx) = fm.shape
#Tmax = times[-1]/tc_D
#L = getGridlen(sdfread(0))
#print(nt,nx)
#ax1.imshow(fm,extent=[0,L,0,Tmax],origin='lower',interpolation='nearest',aspect='auto',cmap='jet')
#xpos=L/20 ; ypos=Tmax/20
#ax1.set_xlabel(r'$x$'+'  '+'[m]',fontsize=18)
#ax1.set_ylabel(r'$t/\tau_{cD}$',fontsize=18)
#ax1.text(xpos,ypos,'(a)',color='k',fontsize=18)

### large sim, zoom
#plot_to_T = int(Tmin/Tmax * nt)
#times = times[0:plot_to_T]
#plot_to_X = int(Lmin/L * nx)
#fm = fm[:plot_to_T,:plot_to_X]
#ax2.imshow(fm,extent=[0,Lmin,0,Tmin],origin='lower',interpolation='nearest',aspect='auto',cmap='jet')
#xpos=Lmin/10 ; ypos=Tmin/10
#ax2.text(xpos,ypos,'(b)',color='k',fontsize=18)
#ax2.set_ylabel(r'$t/\tau_{cD}$',fontsize=18)
#ax2.set_xlabel(r'$x$'+'  '+'[m]',fontsize=18)
##ax2.tick_params('y',labelbottom=False)
#os.chdir('/storage/space2/phrmsf/')
#fig.savefig('fm_Bz_boxes.png')
#fig.savefig('fm_Bz_boxes.eps')
##plt.show()


### weightings
#sim_lst = ['','traceT_0_00','traceT_D_99_T_01','traceT_D_89_T_11','traceT_0_50']
#for sim in sim_lst:
#	sim_loc = getSimulation('/storage/space2/phrmsf/'+sim)
#	print(sim)
#	d0 = sdfread(0)
#	try:
#		ions = getIonSpecies(d0)
#	except:
#		ions = ['Deuterons','Tritons','Alphas']
#	for ion in ions:
#		try:
#			weight = getQuantity1d(d0,'Particles_Weight_'+ion)
#			print(ion, ' :: ', np.mean(weight))
#			print('Npseudo :: '+str(len(weight)))
#		except:
#			print('No weights.')
#	print('end_')

###############################################################################################################
###############################################################################################################

#sims = ['traceT_D_99_T_01','traceT_D_95_T_05','traceT_D_89_T_11','traceT_D_82_T_18','traceT_0_50']
#Lambda = 2e-3
#Z3 = 2
#Mu = np.array([0.01,0.05,0.11,0.18,0.50])
#m=0
#ind_lst = np.arange(12001)
#species = ['Deuterons','Tritons']
#eD = [1.]
#eT = [0.]
#for sim in sims:
#	sim_loc = getSimulation('/storage/space2/phrmsf/'+sim)
#	Energy = np.zeros((len(species),12001))
#	for i in range(len(species)):
#		Energy[i,:] = read_pkl(species[i]+'_KE')
#	times = read_pkl('times')
#	dt = times/len(ind_lst)
#	nD = getMeanquantity(sdfread(0),'Derived_Number_Density_Deuterons')
#	nT = getMeanquantity(sdfread(0),'Derived_Number_Density_Tritons')
#	print(nD,nT)
#	ED = np.sum(Energy[0,:]*dt/nD)
#	ET = np.sum(Energy[1,:]*dt/nT)
#	E = ED+ET
#	print(E,ED,ET)
#	eRat1 = ED/E #Energy[0,:]/Energy[1,:]
#	eRat2 = ET/E #Energy[0,:]/Energy[1,:]
#	eD.append(eRat1)
#	eT.append(eRat2)
#	plt.scatter(Mu[m],eRat1,color='b')
#	plt.scatter(Mu[m],eRat2,color='r')

##	num = 100*(1-Z3*Lambda)-Mu[m]
##	plt.axhline(num/Mu[m],linestyle='--',color='darkgrey')
#	m+=1
#Mu = np.linspace(0,.5,1000)
#Kappa = 1-Mu-Z3*Lambda
#ne = 1e19
#nD = Kappa*ne
#nT = Mu*ne
#mDmT = const.me_to_mD/const.me_to_mT ; Z2 = 1
#mu_eq = (1-Z3*Lambda)/(Z2+mDmT)

##eD.append(0.5)
##eT.append(0.5)
##mu = [0,0.01,0.05,0.11,0.18,0.5,mu_eq]
##plt.scatter(mu,eD,color='b') ; plt.scatter(mu,eT,color='r')

##def cubicfunc(X,a,b,c,d):
##	return c+(b*(a-X)**d)
##popt,pcov=curve_fit(cubicfunc,mu,eD,maxfev=10000)
##ydata = cubicfunc(Mu,popt[0],popt[1],popt[2],popt[3])
##plt.plot(Mu,ydata)
#plt.axvline(mu_eq,linestyle='--',color='darkgrey')
#plt.show()

###############################################################################################################
###############################################################################################################

#sims = ['traceT_0_00','traceT_D_99_T_01','traceT_D_89_T_11','traceT_0_50']
#print('sim','LDe','LDD','va/c','vperp/va')
#for sim in sims:
#	simloc = getSimulation('/storage/space2/phrmsf/'+sim)
#	d0 = sdfread(0)
#	Emin = 3.5e6 * const.qe
#	theta = 0.22*const.PI
#	mass = getMass('Alphas')
#	vz = np.sqrt(2*Emin/mass)
#	v_perp = vz * np.cos(theta)
#	vA = getAlfvenVel(d0)
#	LDe = getDebyeLength(d0,'Electrons')
#	LDD = getDebyeLength(d0,'Deuterons')
#	print([sim,LDe,LDD,vA/const.c,v_perp/vA])

#	wcd = getCyclotronFreq(sdfread(0),'Deuterons')
#	omegas = wcd*np.linspace(0,25,10000)
#	_,k2,_ = coldplasmadispersion(sdfread(0),'Deuterons','',omegas)
#	PIxx, PIxy, PIyy = Chi0Calc(sdfread(0),vz,k2,omegas,species='Deuterons',theta=89)
#	print(PIxx,PIxy,PIyy)

#Z1 = 1
#Z2 = Z3 = 2
#Lambda = np.arange(1e-6,1e-2,1e-7)
#Mu = np.arange(0,1,0.001)
#L,M = np.meshgrid(Lambda,Mu)
#C3 = 15
##CT = ((C3/L))*((1+1/Z1)+M*(1-Z2/Z1)+L*(1-Z3/Z1))
##plt.imshow(np.log10(CT),interpolation='nearest',origin='lower',aspect='auto',extent=[Mu[0],Mu[-1],Lambda[0],Lambda[-1]])
##plt.colorbar()
##plt.show()
#n0 = 5e19
#T0 = 2*11.6e6
#L = 3.
#dx = np.sqrt(const.kb*const.e0*T0/(n0*const.qe**2)) 
#Nx = L/dx
#CT = ((C3/Lambda))*((1+1/Z1)+0.1*(1-Z2/Z1)+Lambda*(1-Z3/Z1))
#NT = CT*Nx
#print(CT,NT)
#plt.plot(Lambda,np.log10(CT),color='k') 
#plt.plot(Lambda,np.log10(NT),color='b') 
#plt.xlabel(r'$\Lambda$')
#plt.ylabel(r'log(Particles) '+'  '+r'$[10^6]$')
#plt.show()

###############################################################################################################
###############################################################################################################




