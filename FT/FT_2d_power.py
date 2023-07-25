
from func_load import *
from matplotlib.gridspec import GridSpec

'''
    Make a multi-panelled plot with FT_2d and power spectra in a 2x2 (technically 2x4) array 
    where each axes has a 1.5:1 width ratio plot of FT_2d:power. Annotated on each is the 
    tritium concentration. Shared x and y-axis according to which plot they represent (links
    back to ax1 and ax2 for 0% case)
    
    Haven't included the color-bar for the FT2d plots, which would be useful, but just can't be
    implemented in an appealing way (TODO: could just do long bar at the top like before?)
'''


#def outside_ticks(fig):
#	for i, ax in enumerate(fig.axes):
#		ax.tick_params(axis='both',direction='out',top=False,right=False,left=True,bottom=True)

def boutside_ticks(lax):
	for ax in lax:
		ax.tick_params(axis='both',direction='out',top=False,right=False,left=True,bottom=True)
	
def xoutside_ticks(lax):
	for ax in lax:
		ax.tick_params(axis='x',direction='out',top=False,right=False,left=False,bottom=True)

def ignorex(lax):
    for ax in lax:
        ax.tick_params(labelbottom=False)

def ignorey(lax):
    for ax in lax:
        ax.tick_params(labelleft=False)

def getFT2d_and_Power(sim):
	simLoc = getSimulation('/storage/space2/phrmsf/'+sim)
	FT_2d = read_pkl('FT_2d_Magnetic_Field_Bz')
	times = read_pkl('times')
	dt = times[-1]/len(times)
	dx = getdxyz(sdfread(0))
	LDe = getDebyeLength(sdfread(0),'Electrons')

	## freq limits
	dw = 2*0.5*const.PI/times[-1]
	klim = 2*0.5*const.PI/dx
	wlim = 2*0.5*const.PI/dt
	vA = getAlfvenVel(sdfread(0))
	wcyc = getCyclotronFreq(sdfread(0),'Deuterons')
	wnorm = wcyc
	knorm = wcyc/vA#1/LDe
	klim_prime = klim/knorm
	wlim_prime = wlim/wnorm
	kmax = 30; wmax = 25
	(nw,nk) = FT_2d.shape
	print('kmax : {} , wmax :{}'.format(kmax,wmax))

	## chopping and plotting
	FT_2d = FT_2d[:int(nw*wmax/wlim_prime),:int(nk*kmax/klim_prime)]
	
	## power spectra plotting
	power = 10**read_pkl('log10_power')
	omegas = read_pkl('omegas_power')
	dw = (omegas[-1]-omegas[0])/len(omegas)
	psd = power/dw
	return FT_2d, psd, omegas, wcyc




sim_lst = ['traceT_0_00','traceT_D_99_T_01','traceT_D_89_T_11','traceT_0_50']
kmax = 30; wmax = 25

fig = plt.figure(figsize=(10,6))
#fig.suptitle("Spatiotemporal Fourier transforms and power spectra")

## 0% & 11%
gs1 = GridSpec(2, 5, left=0.05, right=0.49, wspace=0., hspace=0.05)# right spacing creates space for gs2 
# 0%
#FT2d
FT2d,power0,omegas,wcyc = getFT2d_and_Power(sim_lst[0])
ax1 = fig.add_subplot(gs1[0, :3])
im0=ax1.imshow(np.log10(FT2d),**kwargs,cmap='magma',extent=[0,kmax,0,wmax],vmin=-2,vmax=6)
#power
ax2 = fig.add_subplot(gs1[0, 3:],sharey=ax1)
ax2.plot(power0,omegas/wcyc)
ax2.set_ylim(0,wmax)
ax2.set_xscale('log')
ax2.text(10**4.5,2.5,r'$0\%$') # data coordinates

# 11%
#FT2d 
FT2d,power,omegas,wcyc = getFT2d_and_Power(sim_lst[1])
ax3 = fig.add_subplot(gs1[1, :3],sharey=ax1)
im11=ax3.imshow(np.log10(FT2d),**kwargs,cmap='magma',extent=[0,kmax,0,wmax],vmin=-2,vmax=6)
#power
ax4 = fig.add_subplot(gs1[1, 3:],sharex=ax2,sharey=ax1)
ax4.plot(power,omegas/wcyc)
ax4.set_ylim(0,wmax)
ax4.set_xscale('log')
ax4.text(10**4.5,2.5,r'$11\%$') # data coordinates

#===================================#

## 1% & 50%
gs2 = GridSpec(2, 5, left=0.51, right=0.98, wspace=0., hspace=0.05)
# 1%
#FT2d
FT2d,power,omegas,wcyc = getFT2d_and_Power(sim_lst[2])
ax5 = fig.add_subplot(gs2[0, :3],sharey=ax1)
im1=ax5.imshow(np.log10(FT2d),**kwargs,cmap='magma',extent=[0,kmax,0,wmax],vmin=-2,vmax=6)
#power
ax6 = fig.add_subplot(gs2[0, 3:],sharex=ax2,sharey=ax1)
ax6.plot(power,omegas/wcyc)
ax6.set_ylim(0,wmax)
ax6.set_xscale('log')
ax6.text(10**4.5,2.5,r'$1\%$') # data coordinates

# 50%
#FT2d
FT2d,power,omegas,wcyc = getFT2d_and_Power(sim_lst[3])
ax7 = fig.add_subplot(gs2[1, :3],sharey=ax1)
im50=ax7.imshow(np.log10(FT2d),**kwargs,cmap='magma',extent=[0,kmax,0,wmax],vmin=-2,vmax=6)
#power
ax8 = fig.add_subplot(gs2[1, 3:],sharex=ax2,sharey=ax1)
ax8.plot(power,omegas/wcyc)
ax8.set_ylim(0,wmax)
ax8.set_xscale('log')
ax8.text(10**4.5,2.5,r'$50\%$') # data coordinates


#===================================#
## formatting
boutside_ticks([ax1,ax3,ax5,ax7])
xoutside_ticks([ax2,ax4,ax6,ax8])
#ax_ignorex = [ax1,ax2,ax4,ax5,ax6,ax8]
#ignorex(ax_ignorex)
ax_ignorey = [ax2,ax4,ax5,ax6,ax7,ax8]
ignorey(ax_ignorey)

vals_power=np.array([10,1e2,1e3,1e4,1e5],dtype=float)
labels_power=[r'$10^1$',r'$10^2$',r'$10^3$',r'$10^4$',r'$10^5$']
print(labels_power)
#ax2.set_xticklabels(vals_power)
ax4.set_xticks(vals_power)
ax4.set_xticklabels(labels_power)
#ax6.set_xticklabels(vals_power)
ax8.set_xticks(vals_power)
ax8.set_xticklabels(labels_power)

ax1.set_xticklabels([])
ax5.set_xticklabels([])
ax3.set_xticklabels([0,5,10,15,20,25,30])
ax7.set_xticklabels([0,5,10,15,20,25,30])

ax1.set_ylabel(r'$\omega/\Omega_D$',**tnrfont)
ax3.set_ylabel(r'$\omega/\Omega_D$',**tnrfont)
ax3.set_xlabel(r'$kv_A/\Omega_D$',**tnrfont)
ax7.set_xlabel(r'$kv_A/\Omega_D$',**tnrfont)
ax4.set_xlabel(r'PSD',**tnrfont)
ax8.set_xlabel(r'PSD',**tnrfont)

os.chdir('/storage/space2/phrmsf/paper/')
fig.savefig('FT_2d_and_power.png',bbox_inches='tight')
