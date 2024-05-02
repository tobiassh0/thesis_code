
if __name__=='__main__':
    from func_load import *
    import lib.batch_load as bl
    import time

    homeloc='/home/space/phrmsf/Documents/EPOCH/5_devel/epoch1d/'
    # two particles, unequal masses, (un)equal NDW
    AB = 'AB'
    AB_NDW = 'AB_NDW'
    # three particles, unequal masses, (un)equal NDW
    ABC = 'ABC'
    ABC_NDW = 'ABC_NDW'
    sims = [AB, AB_NDW, ABC, ABC_NDW] # np.sort([i for i in os.listdir(homeloc) if 'AB' in i])

    # # run batch sim script
    # for sim in sims:
    #     print(sim)
    #     ttime = 0
    #     lsdf = list_sdf(getSimulation(homeloc+sim))[-1] # refresh sdf count
    #     while lsdf != 750:
    #         lsdf = list_sdf(getSimulation(homeloc+sim))[-1] # refresh sdf count
    #         ttime += 5 # add total time
    #         print('Waiting for {} to finish...{:.0f}min'.format(sim, ttime))
    #         time.sleep(5*60) # wait 5 minutes
    #     bl.Simulation(homeloc+sim)

    tcp = 2*const.PI/(const.qe*2.1/getMass('Protons'))
    tce = 2*const.PI/(const.qe*2.1/getMass('Electrons'))
    print(tcp,tce,tcp/tce)
    sys.exit()
    fig,axs=plt.subplots(nrows=2,ncols=2,figsize=(10,5),sharey='row',sharex=True)#,layout='constrained')
    fig.subplots_adjust(hspace=0.1,wspace=0.1)
    axs=axs.ravel() # ravel axes
    for i in range(len(sims)):
        simloc=getSimulation(homeloc+sims[i])
        # load times and normalise
        times = read_pkl('times')/tcp
        # load particles and energies
        Bzenergy = read_pkl('Bzenergy')*const.mult_magn_field
        Exenergy = read_pkl('Exenergy')*const.mult_elec_field
        # plot all traces with labels
        axs[i].plot(times,Bzenergy-np.mean(Bzenergy[:10]))
        axs[i].plot(times,Exenergy-np.mean(Exenergy[:10]))
        species = list(filter(None,getIonSpecies(sdfread(0))))
        for spec in species:
            spec_energy = read_pkl(spec+'_KEdens')
            mean_spec_energy = np.mean(spec_energy[:10])
            axs[i].plot(times,spec_energy-mean_spec_energy,label=spec)
        # # annotate each panel
        # axs[i].annotate(sims[i],xy=(0.5,0.5),xycoords='axes fraction',ha='center',va='center')
    # formatting and save
    legend = axs[-1].legend([r'$E_x \epsilon/2$',r'$\Delta B_z^2/2\mu_0$',r'$D_2$',r'$p$',r'$\alpha$'],loc='upper center',\
                        ncol=5,bbox_to_anchor=(-0.05,2.4),borderpad=0.1)
    fig.supylabel(r'$\Delta u \; [Jm^{-3}]$',x=0.04,**tnrfont)
    fig.supxlabel(r'$t/\tau_{cp}$',y=-0.03,**tnrfont)
    axs[0].set_xlim(0,1)
    axs[0].set_ylim(-2,12)
    axs[2].set_ylim(-1,3)
    fig.savefig(homeloc+'two_three_comp.png',bbox_inches='tight')
    # plt.show()


