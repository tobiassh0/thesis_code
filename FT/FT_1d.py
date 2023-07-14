import numpy as np
import matplotlib.pyplot as plt 
import my_constants as const
from list_new import *
from scipy.interpolate import interp2d

sim_lst = ['traceT_0_00','traceT_D_99_T_01','traceT_D_89_T_11','traceT_0_50']
FTs=[]
PowAllk=[]
PowAllt=[]
#fig, axs = plt.subplots(nrows=3,figsize=(12,6),sharex=True) #(8,4)
fig, axs = plt.subplots(nrows=len(sim_lst),figsize=(12,8),sharex=True) #(8,4)
fig.subplots_adjust(hspace=.1)
labels=[r'$0\%$',r'$1\%$',r'$11\%$',r'$50\%$']
for i in range(len(sim_lst)):
	## setup
	print(sim_lst[i])
	simLoc = getSimulation('/storage/space2/phrmsf/'+sim_lst[i])
	FT_1d = read_pkl('FT_1d_Magnetic_Field_Bz')
	times = read_pkl('times')
	dt = times[-1]/len(times)
	dx = getdxyz(sdfread(0))
	LDe = getDebyeLength(sdfread(0),'Electrons')
	## freq limits
	klim = 2*0.5*const.PI/dx
	tlim = times[-1]
	vA = getAlfvenVel(sdfread(0))
	wcyc = getCyclotronFreq(sdfread(0),'Deuterons')
	tnorm = 2*const.PI/wcyc
	knorm = wcyc/vA #1/LDe
	klim_prime = klim/knorm
	tlim_prime = tlim/tnorm
	(nt,nk) = FT_1d.shape
	kmax = 40 #0.035
	tplot = times[-1]/tnorm
	print(nt,nk)
	## chopping and plotting
	FT_1d = FT_1d[1:int(nt*tplot/tlim_prime),:int(nk*kmax/klim_prime)]
	FTs.append(FT_1d)
#	Pow = np.zeros(FT_1d.shape[0])
#	kstart, kstop = 0, int(nk*kmax/klim_prime) #used in case the matrix is not chopped before-hand
#	for i in range(1,FT_1d.shape[0]): # ignore first step through time (will be very large)
#		Pow[i] = np.sum(FT_1d[i,kstart:kstop]**2)
#	PowAllk.append(Pow)
#	Pow = np.zeros(FT_1d.shape[1])
#	for j in range(FT_1d.shape[1]):
#		Pow[j] = np.sum(FT_1d[:,j]**2)
#	PowAllt.append(Pow)
	im = axs[i].imshow(np.log10(FT_1d),interpolation='nearest',cmap='seismic',origin='lower',aspect='auto',extent=[0,kmax,0,tplot],vmin=-3,vmax=4)
	axs[i].set_ylabel(r'$t/\tau_{cD}$',fontsize=22)
	axs[i].annotate(labels[i],xy=[0.9,0.75],xycoords='axes fraction',fontsize=20)#,bbox=dict(fc="white"))
	axs[i].tick_params(axis='both',direction='out',top=False,right=False,left=True,bottom=True)

#	fig.colorbar(im)
#	plt.show()

#axs[len(axs)-1].set_xlabel(r'$k\lambda_{De}$',fontsize=22)
axs[len(axs)-1].set_xlabel(r'$kv_A/\Omega_D$',fontsize=22)#\lambda_{De}

## colorbar
p00 = axs[0].get_position().get_points().flatten()
p01 = axs[1].get_position().get_points().flatten() 
p02 = axs[2].get_position().get_points().flatten()
p03 = axs[3].get_position().get_points().flatten()
ax0_cbar = fig.add_axes([p00[3]+0.01, p03[1], 0.01, p00[3]-p03[1]]) # [left bottom width height]
plt.colorbar(im, cax=ax0_cbar, orientation='vertical')
#fig.text(p00[3]+0.05, (p02[3]), r'$\log(B_z)$', va='center', rotation='vertical',fontsize=20) # xlabel, top row

## save
os.chdir('/storage/space2/phrmsf/paper/remake/')
fig.savefig('FT_1d_Bz_collect_w_50.png')
fig.savefig('FT_1d_Bz_collect_w_50.eps')
#plt.show()









