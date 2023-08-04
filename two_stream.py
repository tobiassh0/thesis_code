from func_load import *
home = os.getcwd()
twostreamloc = '/home/space/phrmsf/Documents/EPOCH/epoch-4.17.16/epoch1d/old/old_sims/twoStream'
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
E_left=np.zeros((len(times),len(v0))) ; E_right=np.zeros((len(times),len(v0)))
E_lhist=np.zeros((len(times),bins)) ; E_rhist=np.zeros((len(times),bins))
Erange = (0,10*Emax)

#fig,ax=plt.subplots(nrows=len(times),figsize=(8,int(1.5*len(times)))) ## uncomment for a nrow plot of histograms
for t, c in zip(times,np.arange(0,len(times),1)):
#	print(t,c)
	# load particle velocities
	vx_left = getQuantity1d(sdfread(t),'Particles_Vx_Left_Electrons')
	vx_right= getQuantity1d(sdfread(t),'Particles_Vx_Right_Electrons')
	E_left[t,:] = 0.5*const.me*vx_left**2
	E_right[t,:]= 0.5*const.me*vx_right**2
	elh, lbin_diff = np.histogram(E_left[t,:],bins=bins,range=Erange,density=True)
	erh, rbin_diff = np.histogram(E_right[t,:],bins=bins,range=Erange,density=True)	
	E_lhist[t,:] = elh*np.diff(lbin_diff)
	E_rhist[t,:] = erh*np.diff(rbin_diff)
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
iml = ax[0].imshow(np.log10(E_lhist),extent=[Erange[0],Erange[1]/(1000*const.qe),0,T],**kwargs,clim=(-6,-2))
imr = ax[1].imshow(np.log10(E_rhist),extent=[Erange[0],Erange[1]/(1000*const.qe),0,T],**kwargs,clim=(-6,-2))
# outside ticks
boutside_ticks(ax)
# colorbar(s) & formatting
p1 = ax[1].get_position().get_points().flatten()
ax1_cbar = fig.add_axes([p1[2]+0.025, p1[1], 0.01, p1[3]-p1[1]]) # [left bottom width height]
rcbar = plt.colorbar(imr, cax=ax1_cbar, orientation='vertical')
ax[0].set_ylabel('Time,  '+r'$\tau_{pe}$',**tnrfont)
labels=['Left Electrons','Right Electrons']
for i in range(len(ax)):
	ax[i].set_xlabel('Energy,  '+r'$keV$',**tnrfont)
	ax[i].set_xticks([0,0.025,0.05,0.075,0.1,0.15,0.2])
	ax[i].set_title(labels[i],**tnrfont)
	ax[i].set_xlim(0,0.12)
print('Normalised such that the sum of histogram in a given time bin (Y-axis) equals 1.')
plt.savefig(home+'/twoStream_Energy_hist.png')
plt.show()

## Fourier transform (in freq) to find freq of oscillation
#plt.clf()
#FT1d = get1dTransform(E_lhist.T)
#print(FT1d.shape)
#dt = tend/len(times)
#wlim = 0.5*2*const.PI/dt
#wpe = getPlasmaFreq(d0,'Left_Electrons')
#plt.imshow(np.log10(FT1d),extent=[0,FT1d.shape[1],0,wlim/wpe],**kwargs) ; plt.show()
