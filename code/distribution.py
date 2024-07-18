
from func_load import *
from scipy import stats
# from sklearn.neighbors import KernelDensity

# load sim
home = '/storage/space2/phrmsf/traceT/'
sim = getSimulation(home+'traceT_D_100_T_00_v3')
ind_lst = list_sdf(sim)
d0 = sdfread(0)
dlast = sdfread(ind_lst[-1])
species = 'Deuterons'

temp = (2/3)*getMeanquantity(d0,'Derived_Average_Particle_Energy_'+species)/const.kb
E0 = (3/2)*const.kb*temp # 3.5e6 * const.qe
mass = getMass(species)
v0 = np.sqrt(2*E0/mass) # vth for bulk ions
vA = getAlfvenVel(d0)
vrange = (0.0,0.1)
evalpoints = np.linspace(vrange[0],vrange[1],1000)
dv = (evalpoints[-1]-evalpoints[0])/len(evalpoints)

fig,ax=plt.subplots(figsize=(7,3),layout='constrained')
darr = [d0,dlast]
colors = ['b','r']
labels=[r'$t/\tau_{ci}=0.0$',r'$t/\tau_{ci}=7.0$']
# loop through time-steps
for i in range(len(darr)):
    # load velocities (equiv. of energies)
    vx = getQuantity1d(darr[i],'Particles_Px_'+species)/mass
    vy = getQuantity1d(darr[i],'Particles_Py_'+species)/mass
    vz = getQuantity1d(darr[i],'Particles_Pz_'+species)/mass
    vi = np.sqrt(vx**2 + vy**2 + vz**2)
    print('loaded velocities')
    # plotting
    ax.hist(vi/vA,density=True,range=vrange,bins=1000,edgecolor='none',alpha=0.5,facecolor=colors[i])
    # # plot dist KDE at start
    # kde= stats.gaussian_kde(vi/vA,bw_method='silverman')
    # yval = kde.pdf(evalpoints)
    # ax.fill_between(evalpoints,yval,0,color=colors[i])
    print('plotted')
    # # check integrals --> 1
    # print('integral ',np.sum(yval*dv))

ax.legend(labels,loc='best',frameon=False)
ax.set_xlabel(r'$|\mathbf{v}|/v_A$',**tnrfont)
ax.set_ylabel(r'$f_i(|\mathbf{v}|)$',**tnrfont)
ax.set_xlim(vrange)
fig.savefig('/home/space/phrmsf/Documents/thesis_code/fi_dist.png')
# plt.show()
