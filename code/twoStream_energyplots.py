
from func_load import *

home = os.getcwd()
#twostreamloc = '/home/space/phrmsf/Documents/EPOCH/epoch-4.17.16/epoch1d/old/old_sims/twoStream_long'
twostreamloc = '/home/space/phrmsf/Documents/EPOCH/5_devel/epoch1d/2stream'
simLoc = getSimulation(twostreamloc)

# files and times
d0 = sdfread(0)
ind_lst = list_sdf(simLoc)
#times = ind_lst[::len(ind_lst)//4]
times = ind_lst

# physical parameters
LDe = getDebyeLength(d0,'Left_Electrons')
dend = sdfread(times[-1])
tau_pe = 2*const.PI/getPlasmaFreq(d0,'Left_Electrons')
tend = dend.__dict__['Header']['time']
v0 = getQuantity1d(d0,'Particles_Vx_Right_Electrons')

# histograms, bins and setup
bins = 5000
Emin, Emax = (np.min(0.5*const.me*v0**2), np.max(0.5*const.me*v0**2))
E_left=np.zeros((len(times),len(v0))) ; E_right=np.zeros((len(times),len(v0))) ; E_ions=np.zeros((len(times),2*len(v0))) # twice as many protons
E_lhist=np.zeros((len(times),bins)) ; E_rhist=np.zeros((len(times),bins)) ; E_ihist=np.zeros((len(times),bins))
Erange = (0,10*Emax)
vx_left = np.zeros((len(times),len(v0))) ; vx_right = np.zeros((len(times),len(v0))) ; vx_ions = np.zeros((len(times),2*len(v0)))

# ------------------------------------------------------- #
## Energy matrices through time (cigarette plots)
#fig,ax=plt.subplots(nrows=len(times),figsize=(8,int(1.5*len(times)))) ## uncomment for a nrow plot of histograms
for t, c in zip(times,np.arange(0,len(times),1)):
#	print(t,c)
	# load particle velocities
	vx_left[t,:]  = getQuantity1d(sdfread(t),'Particles_Vx_Left_Electrons')
	vx_right[t,:] = getQuantity1d(sdfread(t),'Particles_Vx_Right_Electrons')
	vx_ions[t,:]  = getQuantity1d(sdfread(t),'Particles_Vx_Protons')
	# convert to energy
	E_left[t,:] = 0.5*const.me*vx_left[t,:]**2
	E_right[t,:]= 0.5*const.me*vx_right[t,:]**2
	E_ions[t,:] = 0.5*const.mp*vx_ions[t,:]**2
	# calculate histograms
	elh, lbin_diff = np.histogram(E_left[t,:],bins=bins,range=Erange,density=True)
	erh, rbin_diff = np.histogram(E_right[t,:],bins=bins,range=Erange,density=True)	
	eih, ibin_diff = np.histogram(E_ions[t,:],bins=bins,range=Erange,density=True)	
	# append to matrices
	E_lhist[t,:] = elh*np.diff(lbin_diff)
	E_rhist[t,:] = erh*np.diff(rbin_diff)
	E_ihist[t,:] = eih*np.diff(ibin_diff)
	# plot hist of energy at times
#	ax[c].hist(E_left/(1000*const.qe),color='r',bins=100,density=True) ## uncomment for a nrow plot of histograms
#	ax[c].hist(E_right/(1000*const.qe),color='b',bins=100,density=True) ## uncomment for a nrow plot of histograms

# plot energy & hist
fig,ax=plt.subplots(ncols=2,figsize=(8,9),sharey=True)
fig.subplots_adjust(wspace=0.05)
L, T = (getGridlen(d0)/LDe, tend/tau_pe)
# energy
#ax[0].imshow(np.log10(E_left),extent=[Erange[0],Erange[1]/(1000*const.qe),0,T],**kwargs)
#ax[1].imshow(np.log10(E_right),extent=[Erange[0],Erange[1]/(1000*const.qe),0,T],**kwargs)
# histogram
iml = ax[0].imshow(np.log10(E_lhist),extent=[Erange[0],Erange[1]/(1000*const.qe),0,T],**kwargs,clim=(-6,-2),cmap='viridis')
imr = ax[1].imshow(np.log10(E_rhist),extent=[Erange[0],Erange[1]/(1000*const.qe),0,T],**kwargs,clim=(-6,-2),cmap='viridis')
#imi = ax[2].imshow(np.log10(E_ihist),extent=[Erange[0],Erange[1]/(1000*const.qe),0,T],**kwargs,clim=(-6,-2))
# outside ticks
boutside_ticks(ax)
# colorbar(s) & formatting
p1 = ax[1].get_position().get_points().flatten()
ax1_cbar = fig.add_axes([p1[2]+0.025, p1[1], 0.01, p1[3]-p1[1]]) # [left bottom width height]
rcbar = plt.colorbar(imr, cax=ax1_cbar, orientation='vertical')
ax[0].set_ylabel('Time,  '+r'$\tau_{pe}$',**tnrfont)
labels=['Left Electrons','Right Electrons','Protons']
for i in range(len(ax)):
	ax[i].set_xlabel('Energy,  '+r'$keV$',**tnrfont)
	ax[i].set_xticks([0,0.025,0.05,0.075,0.1,0.15,0.2])
	ax[i].set_title(labels[i],**tnrfont)
	ax[i].set_xlim(0,0.12)
print('Normalised such that the sum of histogram in a given time bin (Y-axis) equals 1.')
fig.savefig(home+'/twoStream_Energy_hist.png')
#plt.show()
plt.clf()

# ------------------------------------------------------- #
## Electrostatic potential plots
#phi_Ex = np.zeros((len(ind_lst),))
#for i in ind_lst:

#	phi_Ex = 

# ------------------------------------------------------- #
## Phase space plots
fig,axs=plt.subplots(nrows=4,ncols=1,figsize=(5,10),sharex=True,sharey=True)
fig.subplots_adjust(hspace=0.1,wspace=0.075)
ax=axs.ravel()

plasma_freq = getPlasmaFreq(sdfread(0),'Left_Electrons')
plasma_prd  = 2*const.PI/plasma_freq
ttimes = np.array([0,2,4,6])
ftimes = [int(len(times)*(ttimes[i]/(tend/plasma_prd))) for i in range(len(ttimes))]
ftimes[-1] = 100
frac = 1
temp = 273 ; dens = 10
v_the = np.sqrt(2*const.kb*temp/const.me)
color = ['r','b','g']
species = ['Left_Electrons','Right_Electrons'] #getAllSpecies(sdfread(0))

for i in range(len(ax)):
	label = r'$t\omega_{pe}/2\pi=$'+str(np.around(ttimes[i],1))
	for j in range(len(species)):
		vx_name = sdfread(ftimes[i]).__dict__['Particles_Vx_'+species[j]]
		x_name = vx_name.grid.data[0]/LDe 		# Normalises distance to units of Debye length
		ax[i].scatter(x_name,vx_name.data/v_the,color=color[j],s=1.) 		# Normalises speed to v_thermal at t=0
	ax[i].text(100.,90.,label,color='black',fontsize=16,bbox=dict(facecolor='white',edgecolor='black',boxstyle='square,pad=0.25'))
	ax[i].set_xlim(0,1400)
	ax[i].set_ylim(-100,100)
	ax[i].set_ylabel(r'$v_x/v_{e,th}$',fontsize=20)
ax[-1].set_xticks([0,350,700,1050,1400])
ax[-1].set_xlabel(r'$x/\lambda_{De}$',fontsize=20)
print(home)
fig.savefig(home+'/twoStream_phase.png',bbox_inches='tight')
#plt.show()

