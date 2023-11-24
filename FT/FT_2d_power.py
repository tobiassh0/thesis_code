
from func_load import *
from matplotlib.gridspec import GridSpec
import FT_2d_vMCI as vmc

'''
    Make a multi-panelled plot with FT_2d and power spectra in a 2x2 (technically 2x4) array 
    where each axes has a 1.5:1 width ratio plot of FT_2d:power. Annotated on each is the 
    tritium concentration. Shared x and y-axis according to which plot they represent (links
    back to ax1 and ax2 for 0% case)
    
    Haven't included the color-bar for the FT2d plots, which would be useful, but just can't be
    implemented in an appealing way (TODO: could just do long bar at the top like before?)
'''


def getFT2d_and_Power(sim,kmax=30,wmax=20):
	simLoc = getSimulation('/storage/space2/phrmsf/traceT/'+sim)
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



width_FT2d = 3

sim_lst = ['traceT_0_00','traceT_D_99_T_01','traceT_D_89_T_11','traceT_0_50']
kmax = 50; wmax = 30
fig = plt.figure(figsize=(11,6))
#fig.suptitle("Spatiotemporal Fourier transforms and power spectra")

# phase and group velocities of strongest MCI wave-packet
xi2 = [0.0,0.01,0.11,0.50]
v_ph, v_gr, w_MCI, k_MCI = vmc.v_MCI(['/storage/space2/phrmsf/traceT/'+i for i in sim_lst], xi2, plot=False)
vA = []
for i in sim_lst:
	getSimulation('/storage/space2/phrmsf/traceT/'+i)
	vA.append(getAlfvenVel(sdfread(0)))
vA = np.array(vA)
dk = 15
#print(v_gr/const.c)
c = w_MCI-v_gr*k_MCI
k2 = k_MCI + dk ; k1 = k_MCI - dk
w2 = v_gr*k2 + c ; w1 = v_gr*k1 + c

## 0% & 11%
gs1 = GridSpec(2, 5, left=0.05, right=0.49, wspace=0., hspace=0.05)# right spacing creates space for gs2 
# 0%
#FT2d
FT2d,power0,omegas,wcyc = getFT2d_and_Power(sim_lst[0],kmax=kmax,wmax=wmax)
ax1 = fig.add_subplot(gs1[0, :width_FT2d])
im0=ax1.imshow(np.log10(FT2d),**kwargs,cmap='magma',extent=[0,kmax,0,wmax],vmin=-2,vmax=6)
ax1.plot([k1[0],k2[0]],[w1[0],w2[0]],color='k',linestyle='--')
ax1.text(k_MCI[0],w_MCI[0]-2,r'$v_{gr}/v_A=$'+'{}'.format(np.around(v_gr[0],2)),ha='left',va='top',**tnrfont)
#power
ax2 = fig.add_subplot(gs1[0, width_FT2d:],sharey=ax1)
ax2.plot(power0,omegas/wcyc)
ax2.set_ylim(0,wmax)
ax2.set_xscale('log')
ax2.text(10**4.2,2.5,r'$0\%$',**tnrfont) # data coordinates

# 11%
#FT2d 
FT2d,power,omegas,wcyc = getFT2d_and_Power(sim_lst[2],kmax=kmax,wmax=wmax)
ax3 = fig.add_subplot(gs1[1, :width_FT2d],sharey=ax1)
im11=ax3.imshow(np.log10(FT2d),**kwargs,cmap='magma',extent=[0,kmax,0,wmax],vmin=-2,vmax=6)
ax3.plot([k1[2],k2[2]],[w1[2],w2[2]],color='k',linestyle='--')
ax3.text(k_MCI[2],w_MCI[2]-2,r'$v_{gr}/v_A=$'+'{}'.format(np.around(v_gr[2],2)),ha='left',va='top',**tnrfont)
#power
ax4 = fig.add_subplot(gs1[1, width_FT2d:],sharex=ax2,sharey=ax1)
ax4.plot(power,omegas/wcyc)
ax4.set_ylim(0,wmax)
ax4.set_xscale('log')
ax4.text(10**4.2,2.5,r'$11\%$',**tnrfont) # data coordinates

#===================================#

## 1% & 50%
gs2 = GridSpec(2, 5, left=0.51, right=0.98, wspace=0., hspace=0.05)
# 1%
#FT2d
FT2d,power,omegas,wcyc = getFT2d_and_Power(sim_lst[1],kmax=kmax,wmax=wmax)
ax5 = fig.add_subplot(gs2[0, :width_FT2d],sharey=ax1)
im1=ax5.imshow(np.log10(FT2d),**kwargs,cmap='magma',extent=[0,kmax,0,wmax],vmin=-2,vmax=6)
ax5.plot([k1[1],k2[1]],[w1[1],w2[1]],color='k',linestyle='--')
ax5.text(k_MCI[1],w_MCI[1]-2,r'$v_{gr}/v_A=$'+'{}'.format(np.around(v_gr[1],2)),ha='left',va='top',**tnrfont)
#power
ax6 = fig.add_subplot(gs2[0, width_FT2d:],sharex=ax2,sharey=ax1)
ax6.plot(power,omegas/wcyc)
ax6.set_ylim(0,wmax)
ax6.set_xscale('log')
ax6.text(10**4.2,2.5,r'$1\%$',**tnrfont) # data coordinates

# 50%
#FT2d
FT2d,power,omegas,wcyc = getFT2d_and_Power(sim_lst[3],kmax=kmax,wmax=wmax)
ax7 = fig.add_subplot(gs2[1, :width_FT2d],sharey=ax1)
im50=ax7.imshow(np.log10(FT2d),**kwargs,cmap='magma',extent=[0,kmax,0,wmax],vmin=-2,vmax=6)
ax7.plot([k1[3],k2[3]],[w1[3],w2[3]],color='k',linestyle='--')
ax7.text(k_MCI[3],w_MCI[3]-2,r'$v_{gr}/v_A=$'+'{}'.format(np.around(v_gr[3],2)),ha='left',va='top',**tnrfont)
#power
ax8 = fig.add_subplot(gs2[1, width_FT2d:],sharex=ax2,sharey=ax1)
ax8.plot(power,omegas/wcyc)
ax8.set_ylim(0,wmax)
ax8.set_xscale('log')
ax8.text(10**4.2,2.5,r'$50\%$',**tnrfont) # data coordinates


#===================================#
## formatting
boutside_ticks([ax1,ax3,ax5,ax7])
xoutside_ticks([ax2,ax4,ax6,ax8])
#ax_ignorex = [ax1,ax2,ax4,ax5,ax6,ax8]
#ignorex(ax_ignorex)
ax_ignorey = [ax2,ax4,ax5,ax6,ax7,ax8]
ignorey(ax_ignorey)

# power tick labels
vals_power=np.array([10,1e2,1e3,1e4,1e5],dtype=float)
labels_power=[r'$10^1$',r'$10^2$',r'$10^3$',r'$10^4$',r'$10^5$']
print(labels_power)
#ax2.set_xticklabels(vals_power)
ax4.set_xticks(vals_power)
ax4.set_xticklabels(labels_power,fontsize=14)
#ax6.set_xticklabels(vals_power)
ax8.set_xticks(vals_power)
ax8.set_xticklabels(labels_power,fontsize=14)

# tick labels
ax1.set_xticklabels([])
ax5.set_xticklabels([])
karr = np.arange(0,kmax+10,10)
ax3.set_xticklabels(karr)
ax7.set_xticklabels(karr)

# axis labels
ax1.set_ylabel(r'$\omega/\Omega_D$',**tnrfont)
ax3.set_ylabel(r'$\omega/\Omega_D$',**tnrfont)
ax3.set_xlabel(r'$kv_A/\Omega_D$',**tnrfont)
ax7.set_xlabel(r'$kv_A/\Omega_D$',**tnrfont)
ax4.set_xlabel(r'PSD',**tnrfont)
ax8.set_xlabel(r'PSD',**tnrfont)

# color bar
p0 = ax1.get_position().get_points().flatten()
p6 = ax6.get_position().get_points().flatten()
ax0_cbar = fig.add_axes([p0[0], 0.97, p6[2]-p0[0], 0.02]) # [left bottom width height]
plt.colorbar(im0, cax=ax0_cbar, orientation='horizontal')

plt.show()
os.chdir('/storage/space2/phrmsf/traceT/paper/')
fig.savefig('FT_2d_and_power.png',bbox_inches='tight')
