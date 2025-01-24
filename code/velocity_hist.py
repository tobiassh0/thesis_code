
from func_load import *
import energy.LarmorRadii as lr
from scipy import stats

def getRowsCols(num):
    N = 0
    M = num
    while N+1<=M:
        N+=1
        if num % N == 0:
            M = num // N
    return N, M

def root_find(x,y,search_between=(6.45,6.5)):
    minsearch, maxsearch = search_between
    thresh = (minsearch < x) & (x < maxsearch)
    ythresh = y[thresh]
    yminval = ythresh[np.argmin(ythresh**2)]
    minyarg = np.where(y==yminval)[0][0]
    return x[minyarg], y[minyarg]

def plot_velocity_startend(home,sim,species,vrange=(0.8,1.8),colors=['b','r']):
    ind_lst = list_sdf(sim)
    # load start and end
    d0 = sdfread(0)
    dlast = sdfread(ind_lst[-1])
    darr = [d0,dlast]
    # mass [kg]
    mass = getMass(species)
    # if thermal ions
    if species in ['Deuterons','Tritons','He3']: 
        temp = (2/3)*getMeanquantity(d0,'Derived_Average_Particle_Energy_'+species)/const.kb
        E0 = 3.5e6 * const.qe # (3/2)*const.kb*temp 
        v0 = np.sqrt(2*E0/mass) # vth for bulk ions
        vA = getAlfvenVel(d0)

    # setup fig
    fig,ax=plt.subplots(figsize=(7,3),layout='constrained')

    # loop through time-steps
    for i in range(len(darr)):
        # load velocities (equiv. of energies)
        vx = getQuantity1d(darr[i],'Particles_Px_'+species)/mass
        vy = getQuantity1d(darr[i],'Particles_Py_'+species)/mass
        vz = getQuantity1d(darr[i],'Particles_Pz_'+species)/mass
        vi = np.sqrt(vx**2 + vy**2 + vz**2)
        print('loaded velocities')
        # plot dist KDE at start
        kde= stats.gaussian_kde(vi/vA,bw_method='silverman')
        yval = kde.pdf(evalpoints)
        ax.fill_between(evalpoints,yval,0,color=colors[i],alpha=0.5)
        # # plotting
        # ax.hist(vi/vA,density=True,range=vrange,bins=1000,edgecolor='none',alpha=0.5,facecolor=colors[i])
        # # check integrals --> 1
        # print('integral ',np.sum(yval*(vrange[-1]-vrange[0])/(len(vrange))*vA)))

    # format fig
    ax.legend(labels,loc='best',frameon=False)
    ax.set_xlabel(r'$|\mathbf{v}|/v_A$',**tnrfont)
    ax.set_ylabel(r'$f_\alpha(|\mathbf{v}|)$',**tnrfont)
    ax.set_xlim(vrange)
    fig.savefig(home+'fv_dist_{}.png'.format(species))
    plt.show()

def plotVelocityKDE(home,sim,species=['Protons'],colors=['r'],dims=['x','y','z'],name='',logname='',num_bins=200):
    simloc = getSimulation(home+sim)
    times = read_pkl('times')
    restart_files = lr.para_check_restart(simloc)
    # load first file and parameters
    d0 = sdfread(0)
    B0 = getMeanField3D(d0,'Magnetic_Field_B')
    theta,_ = getMagneticAngle(d0) # radians
    vA = getAlfvenVel(d0)
    tcspec = 2*const.PI/getCyclotronFreq(d0,species[-1])
    # setup one fig
    Nrows,Ncols=getRowsCols(len(restart_files))
    fig,ax=plt.subplots(figsize=(3*Nrows,3*Ncols),nrows=Nrows,ncols=Ncols,sharey=True,sharex=True,layout='constrained')
    ax=ax.ravel()
    # root v solutions
    rootv = np.zeros(len(restart_files))
    # loop through each time
    for i in range(len(restart_files)):
        V=[]
        print(i,restart_files[i])
        # loop through species:
        for k in range(len(species)):
            # loop through dimensions
            for j in range(len(dims)):
                rstfl = sdfread(restart_files[i])
                V.append(getQuantity1d(rstfl,'Particles_P'+dims[j]+'_'+species[k])/getMass(species[k]))
        V=np.array(V).flatten()
        ax[i].hist(V/vA,bins=num_bins,density=True,range=(-1.5,1.5),edgecolor='none')
        # # KDE on total velocities
        # Vkde = stats.gaussian_kde(V/vA)
        # varr = np.linspace(-2*vA,2*vA,num_bins) # np.min(V),np.max(V)
        # dist = Vkde(varr/vA)
        # ax[i].plot(varr/vA,dist,color=colors[k])
        """
        # get gradient of dist
        grad_dist = np.gradient(dist,(varr[-1]-varr[0])/(vA*num_bins))
        ax[i].plot(varr/vA,grad_dist,color='k',zorder=2)
        # minimisation of gradient
        vmin, grad_min = root_find(varr/vA,grad_dist)
        rootv[i] = vmin*vA
        ax[i].scatter(vmin,grad_min,color='r',edgecolor='none',zorder=3)
        """
        # gridlines
        ax[i].grid(True)
        # ax[i].axhline(0,linestyle='--',color='darkgrey',zorder=1)
        ax[i].annotate(r'${:.2f}$'.format(times[restart_files[i]]/tcspec)+r'$\tau_{c\sigma}$',xy=(0.95,0.95),\
                        xycoords='axes fraction',ha='right',va='top',**tnrfont)
    # plot save
    figname = '_'.join([str(len(species)),name,logname,'.pdf'])
    fig.supxlabel(r'$\mathbf{v}/v_A$',**tnrfont)
    fig.supylabel('KDE',**tnrfont) # (green) ; Gradient (black)
    # ax[-1].set_xlim(6.4,6.6)
    fig.savefig(home+'V_KDEgrad_'+figname)
    plt.clf()
    # plt.plot(restart_files,rootv,'-o') ; plt.show()
    return None

def plotVelocitydimKDE(home,sim,species='Protons',color='r',dims=['x','y','z'],name='',logname='',num_bins=200):
    simloc = getSimulation(home+sim)
    times = read_pkl('times')
    restart_files = lr.para_check_restart(simloc)
    # load first file and parameters
    d0 = sdfread(0)
    B0 = getMeanField3D(d0,'Magnetic_Field_B')
    theta,_ = getMagneticAngle(d0) # radians
    vA = getAlfvenVel(d0)
    tcspec = 2*const.PI/getCyclotronFreq(d0,species)
    # loop figs for each dimension
    figarr=[] 
    axarr=[] # shape [dims, restart_files]
    Nrows,Ncols=getRowsCols(len(restart_files))
    for dim in dims:
        # open one fig
        figi,axi=plt.subplots(figsize=(3*Nrows,3*Ncols),nrows=Nrows,ncols=Ncols,sharey=True,sharex=True,layout='constrained')
        axi=axi.ravel()
        # append to array of figs
        figarr.append(figi)
        axarr.append(axi)
    figarr=np.array(figarr)
    axarr=np.array(axarr)
    # loop through each time
    for i in range(len(restart_files)):
        V2=0
        print(i,restart_files[i])
        # loop through dimensions
        for j in range(len(dims)):
            rstfl = sdfread(restart_files[i])
            axarr[j,i].axhline(0,linestyle='--',color='darkgrey')
            Vdim=getQuantity1d(rstfl,'Particles_P'+dims[j]+'_'+species)/getMass(species)
            # KDEs # TODO; make this quicker
            Vdimkde = stats.gaussian_kde(Vdim/vA)
            varr = np.linspace(min(Vdim),max(Vdim),num_bins)
            dist = Vdimkde(varr/vA)
            axarr[j,i].plot(varr/vA,dist,color=color) # /np.max(dist)
            # get gradient of dist
            grad_dist = np.gradient(dist,(varr[-1]-varr[0])/(vA*num_bins))
            axarr[j,i].plot(varr/vA,grad_dist,color='k') # /np.max(grad_dist)
            axarr[j,i].annotate(r'${:.2f}$'.format(times[restart_files[i]]/tcspec)+r'$\tau_{c\sigma}$',xy=(0.95,0.95),\
                            xycoords='axes fraction',ha='right',va='top',**tnrfont)
    # plot save
    plt.show()
    figname = '_'.join([species,name,logname,'.pdf'])
    for j in range(len(dims)):
        figarr[j].supxlabel(r'$v_{}/v_A$'.format(dims[j]),**tnrfont)
        figarr[j].supylabel('KDE (green) ; Gradient (black)',**tnrfont)
        figarr[j].savefig(home+'V{}_KDEgrad_'.format(dims[j])+figname)
    return None

def plotVelocityHist(home,sims,species,minspec='Protons',binrange=None,colors=['b','r','g'],labels=[''],num_bins=200,logyscale=True,KDE=True):
    # if binrange == None:
    #     binrange = (-0.15,0.15)
    if len(colors)!=len(species):
        colors = plt.cm.rainbow(np.linspace(0,1,len(species)))
    if len(labels)!=len(sims):
        labels=['']*len(sims)

    for j in range(len(sims)):
        print(sims[j]+' :: {:.2f}%'.format(100*j/len(sims))) # index as %
        simloc = getSimulation(home+sims[j])
        times = read_pkl('times')
        restart_files = lr.para_check_restart(simloc)
        # load first file and parameters
        d0 = sdfread(0)
        B0 = getMeanField3D(d0,'Magnetic_Field_B')
        theta,_ = getMagneticAngle(d0) # radians
        vA = getAlfvenVel(d0)
        tcmin = 2*const.PI/getCyclotronFreq(d0,minspec)
        # setup figures
        Nrows,Ncols=getRowsCols(len(restart_files))
        figx,axx=plt.subplots(figsize=(3*Nrows,3*Ncols),nrows=Nrows,ncols=Ncols,sharey=True,sharex=True,layout='constrained')
        figy,axy=plt.subplots(figsize=(3*Nrows,3*Ncols),nrows=Nrows,ncols=Ncols,sharey=True,sharex=True,layout='constrained')
        figz,axz=plt.subplots(figsize=(3*Nrows,3*Ncols),nrows=Nrows,ncols=Ncols,sharey=True,sharex=True,layout='constrained')
        axx = axx.ravel() ; axy = axy.ravel() ; axz = axz.ravel() # collapse into 1d arrays
        # loop through times
        for i in range(len(restart_files)):
            print(i,restart_files[i])
            # loop through species
            for s in range(len(species)):
                rstfl = sdfread(restart_files[i])
                Vx=getQuantity1d(rstfl,'Particles_Px_'+species[s])/getMass(species[s])
                Vy=getQuantity1d(rstfl,'Particles_Py_'+species[s])/getMass(species[s])
                Vz=getQuantity1d(rstfl,'Particles_Pz_'+species[s])/getMass(species[s])
                # histograms
                _=axx[i].hist(Vx/vA,bins=num_bins,color=colors[s],edgecolor='none',alpha=0.35,density=True,range=binrange)
                _=axy[i].hist(Vy/vA,bins=num_bins,color=colors[s],edgecolor='none',alpha=0.35,density=True,range=binrange)
                _=axz[i].hist(Vz/vA,bins=num_bins,facecolor=colors[s],edgecolor='none',alpha=0.35,density=True,range=binrange)
                # plot KDE on-top
                if KDE:
                    Vxkde = stats.gaussian_kde(Vx/vA)
                    Vykde = stats.gaussian_kde(Vy/vA)
                    Vzkde = stats.gaussian_kde(Vz/vA)
                    varrx = np.linspace(min(Vx),max(Vx),num_bins)
                    varry = np.linspace(min(Vy),max(Vy),num_bins)
                    varrz = np.linspace(min(Vz),max(Vz),num_bins)
                    axx[i].plot(varrx/vA,Vxkde(varrx/vA),color=colors[s])
                    axy[i].plot(varry/vA,Vykde(varry/vA),color=colors[s])
                    axz[i].plot(varrz/vA,Vzkde(varrz/vA),color=colors[s])
            # annotate times
            axx[i].annotate(r'${:.2f}$'.format(times[restart_files[i]]/tcmin)+r'$\tau_{cmin}$',xy=(0.95,0.95),xycoords='axes fraction',ha='right',va='top',**tnrfont)
            axy[i].annotate(r'${:.2f}$'.format(times[restart_files[i]]/tcmin)+r'$\tau_{cmin}$',xy=(0.95,0.95),xycoords='axes fraction',ha='right',va='top',**tnrfont)
            axz[i].annotate(r'${:.2f}$'.format(times[restart_files[i]]/tcmin)+r'$\tau_{cmin}$',xy=(0.95,0.95),xycoords='axes fraction',ha='right',va='top',**tnrfont)
        logname = ''
        if logyscale:
            # log yscale
            axx[0].set_yscale('log') # sharey=True so will account for all
            axy[0].set_yscale('log')
            axz[0].set_yscale('log')
            logname='_log'
        # vx plot save
        figx.supxlabel(r'$v_x/v_A$',**tnrfont)
        figx.supylabel('Normalised Count',**tnrfont)
        figx.savefig(home+'Vx_hist_'+labels[j]+logname+'.pdf')
        # vy plot save
        figy.supxlabel(r'$v_y/v_A$',**tnrfont)
        figy.supylabel('Normalised Count',**tnrfont)
        figy.savefig(home+'Vy_hist_'+labels[j]+logname+'.pdf')
        # vz plot save
        figz.supxlabel(r'$v_z/v_A$',**tnrfont)
        figz.supylabel('Normalised Count',**tnrfont)
        figz.savefig(home+'Vz_hist_'+labels[j]+logname+'.pdf')
        plt.close('all')
    return None

def plotVelocityScatter(home,sim,species='Protons'):
    simloc = getSimulation(home+sim)
    times = read_pkl('times')
    restart_files = lr.para_check_restart(simloc)
    mass = getMass(species)
    vA = getAlfvenVel(sdfread(0))
    tcspec = 2*const.PI/getCyclotronFreq(sdfread(0),species)
    az=45 ; el=30 # degrees
    for i in range(len(restart_files)):
        fig=plt.figure()
        ax=fig.add_subplot(projection='3d')
        print('{:.2f}%'.format(100*i/len(restart_files))) # print progress
        d = sdfread(restart_files[i])
        Vx = getQuantity1d(d,'Particles_Px_'+species)/mass
        Vy = getQuantity1d(d,'Particles_Py_'+species)/mass
        Vz = getQuantity1d(d,'Particles_Pz_'+species)/mass
        ax.plot(Vx/vA, Vy/vA, Vz/vA, ',k', alpha=0.25)
        ax.view_init(elev=el,azim=az)
        ax.set_title('{:.2f}'.format(times[restart_files[i]]/tcspec,**tnrfont)+r'$\tau_{c\sigma}$')
        ax.set_xlabel(r'$v_x/v_A$',**tnrfont)
        ax.set_ylabel(r'$v_y/v_A$',**tnrfont)
        ax.set_zlabel(r'$v_z/v_A$',**tnrfont)
        fig.savefig(home+'VelocityScatter_'+sim+'_'+str(int(restart_files[i]))+'.png')
    plt.close('all') # close all open figs
    return None

if __name__=='__main__':
    # DT
    home = '/storage/space2/phrmsf/traceT/'
    sim = getSimulation(home+'traceT_D_100_T_00_v3')
    _labels=[r'$t/\tau_{c\alpha}=0.0$',r'$t/\tau_{c\alpha}=7.0$']
    # plot velocity hist of Alphas at start and end
    plot_velocity_startend(home,sim,'Alphas',labels=_labels)
    sys.exit()

    # DHe3 
    home = '/storage/space2/phrmsf/lowres_D_He3/'
    sims = np.sort([i for i in os.listdir(home) if 'p_90' in i and '.pdf' not in i])[1:] # remove 0%
    xiHe3 = [str(int(i.split('_')[1]))+'%' for i in sims] # 'Protons_'+
    species=['Deuterons','He3','Protons']
    colors=['b','g','r']
    for i in range(len(sims)):
        print(sims[i])
        plotVelocityKDE(home,sims[i],species=species,colors=colors,dims=['x'],name=xiHe3[i]+'_bulkandring',num_bins=1000)
    sys.exit()
    # plotVelocityScatter(home,sims[0],species='Protons')
    # plotVelocitydimKDE(home,sims[0],species='Protons',color='g',dims=['z'],name='5%',num_bins=1000)
    # plotVelocityHist(home,sims,['Protons'],minspec='Protons',colors=['g'],labels=xiHe3,num_bins=200,logyscale=True)
