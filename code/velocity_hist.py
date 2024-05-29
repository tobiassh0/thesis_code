
from func_load import *
import energy.LarmorRadii as lr
from scipy import stats

""" TODO; parallelise loading velocities for each component, species and time
def para_velocity_get(i,species=['Deuterons','He3']):
    rstfl = sdfread(i)
    VARR = []
    for spec in species:
        Vx=getQuantity1d(rstfl,'Particles_Px_'+spec)/getMass(spec)
        Vy=getQuantity1d(rstfl,'Particles_Py_'+spec)/getMass(spec)
        Vz=getQuantity1d(rstfl,'Particles_Pz_'+spec)/getMass(spec)
        VARR.append([Vx,Vy,Vz])
    return

def para_velocity_hist(restart_files):
    # pass indice of restart and both species, return vel component arrays size [nspec, ntimes, ndim, nbins] = [2, 25, 3, 100]
    pool=mp.Pool(mp.cpu_count())
    arr = np.array(pool.map_async(para_velocity_get,restart_files).get(99999))
    pool.close()
    pass
"""

def getRowsCols(num):
    N = 0
    M = num
    while N+1<=M:
        N+=1
        if num % N == 0:
            M = num // N
    return N, M

def plotVelocityHist(home,sims,species,minspec='Protons',binrange=None,colors=['b','r','g'],labels=[''],num_bins=200,logyscale=True):
    # if binrange == None:
    #     binrange = (-0.15,0.15)
    if len(colors)!=len(species):
        colors = plt.cm.rainbow(np.linspace(0,1,len(species)))
    if len(labels)!=len(sims):
        labels=['']*len(sims)

    for j in range(len(sims)):
        print(sim+' :: {:.2f}%'.format(100*j/len(sims))) # index as %
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
                # # KDEs # TODO; make this quicker
                # Vxkde = stats.gaussian_kde(Vx)
                # Vykde = stats.gaussian_kde(Vy)
                # Vzkde = stats.gaussian_kde(Vz)
                # axx[i].plot(varr,Vxkde(varr),color=colors[s])
                # axy[i].plot(varr,Vykde(varr),color=colors[s])
                # axz[i].plot(varr,Vzkde(varr),color=colors[s])
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
    home = '/storage/space2/phrmsf/lowres_D_He3/'
    sims = np.sort([i for i in os.listdir(home) if 'p_90' in i and '.pdf' not in i])[1:] # remove 0%
    xiHe3 = [i.split('_')[1] for i in sims]
    # plotVelocityScatter(home,sims[0],species='Protons')
    plotVelocityHist(home,sims,['Protons'],minspec='Protons',colors=['g'],labels=xiHe3,num_bins=200,logyscale=True)
    sys.exit()

    binrange=(-0.3,0.3)
    N = 400
    varr = np.linspace(binrange[0],binrange[1],N)
    # loop through xiHe3, times, and species, get velocity components and plot as hist
    species = ['Deuterons','He3']
    colors = ['r','b']
    for sim in sims:
        xiHe3 = int(sim.split('_')[1])
        print('SIM :: {}%'.format(xiHe3))
        simloc = getSimulation(home+sim)
        times = read_pkl('times')
        restart_files = lr.para_check_restart(simloc)
        # restart_files = [0]#,3000,6000,9000,12000]
        d0 = sdfread(0)
        B0 = getMeanField3D(d0,'Magnetic_Field_B')
        theta,_ = getMagneticAngle(d0) # radians
        vA = getAlfvenVel(d0)
        tcp = 2*const.PI/getCyclotronFreq(d0,'Protons')
        # setup figures
        figx,axx=plt.subplots(figsize=(15,15),nrows=5,ncols=5,sharey=True,sharex=True,layout='constrained')
        figy,axy=plt.subplots(figsize=(15,15),nrows=5,ncols=5,sharey=True,sharex=True,layout='constrained')
        figz,axz=plt.subplots(figsize=(15,15),nrows=5,ncols=5,sharey=True,sharex=True,layout='constrained')
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
                _=axx[i].hist(Vx/vA,bins=N,color=colors[s],edgecolor='none',alpha=0.35,density=True,range=binrange)
                _=axy[i].hist(Vy/vA,bins=N,color=colors[s],edgecolor='none',alpha=0.35,density=True,range=binrange)
                _=axz[i].hist(Vz/vA,bins=N,facecolor=colors[s],edgecolor='none',alpha=0.35,density=True,range=binrange)
                # # KDEs # TODO; make this quicker
                # Vxkde = stats.gaussian_kde(Vx)
                # Vykde = stats.gaussian_kde(Vy)
                # Vzkde = stats.gaussian_kde(Vz)
                # axx[i].plot(varr,Vxkde(varr),color=colors[s])
                # axy[i].plot(varr,Vykde(varr),color=colors[s])
                # axz[i].plot(varr,Vzkde(varr),color=colors[s])
            # annotate times
            axx[i].annotate(r'${:.2f}$'.format(times[restart_files[i]]/tcp)+r'$\tau_{cp}$',xy=(0.95,0.95),xycoords='axes fraction',ha='right',va='top',**tnrfont)
            axy[i].annotate(r'${:.2f}$'.format(times[restart_files[i]]/tcp)+r'$\tau_{cp}$',xy=(0.95,0.95),xycoords='axes fraction',ha='right',va='top',**tnrfont)
            axz[i].annotate(r'${:.2f}$'.format(times[restart_files[i]]/tcp)+r'$\tau_{cp}$',xy=(0.95,0.95),xycoords='axes fraction',ha='right',va='top',**tnrfont)
            # log yscale
            axx[i].set_yscale('log')
            axy[i].set_yscale('log')
            axz[i].set_yscale('log')
        # vx plot save
        figx.supxlabel(r'$v_x/v_A$',**tnrfont)
        figx.supylabel('Normalised Count',**tnrfont)
        figx.savefig(home+'Vx_hist_'+str(xiHe3)+'%_log.pdf')
        # vy plot save
        figy.supxlabel(r'$v_y/v_A$',**tnrfont)
        figy.supylabel('Normalised Count',**tnrfont)
        figy.savefig(home+'Vy_hist_'+str(xiHe3)+'%_log.pdf')
        # vz plot save
        figz.supxlabel(r'$v_z/v_A$',**tnrfont)
        figz.supylabel('Normalised Count',**tnrfont)
        figz.savefig(home+'Vz_hist_'+str(xiHe3)+'%_log.pdf')
    # plt.show()

