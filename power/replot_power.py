

if __name__=='__main__':
    from func_load import *

    # load sim
    """
    user input
    """
    home = '/storage/space2/phrmsf/p_B11'
    simloc = getSimulation(home)
    d0 = sdfread(0)
    times = read_pkl('times')

    # normalisation and limits
    minspec = 'Alphas'
    wnorm = getCyclotronFreq(d0,minspec)
    vA = getAlfvenVel(d0)
    wmax = 0.5*2*const.PI/getdt(times)
    dw = 2*wmax/len(times) # account for +ve and -ve wmax covered by N_t files
    """
    user input
    """
    wlim = 20
    load = True

    if load: # load power_pkl file
        log10_power = read_pkl('log10_power')
        omegas = read_pkl('omegas_power')
    else: # re-calculate
        kmax = 0.5*2*const.PI/getdxyz(d0)
        knorm = wnorm/vA
        omegas,log10_power=power(wnorm,wklims=[wmax/wnorm,kmax/knorm],wkmax=[40,100],norm_omega=getOmegaLabel(minspec),quantity='Magnetic_Field_Bz',plot=False,read=False,dump=False,outp=True)

    power = 10**log10_power
    thresh = omegas/wnorm < wlim
    fig,ax=plt.subplots(figsize=(8,5))
    # plot harmonics
    for i in range(wlim+1):
        ax.axvline(i,color='darkgrey',linestyle='--')
    # plot power
    ax.plot(omegas[thresh]/wnorm,power[thresh]/dw) # freq, psd
    # formatting
    ax.set_ylabel('PSD',**tnrfont)
    ax.set_xlabel(r'$\omega/$'+getOmegaLabel(minspec),**tnrfont)
    ax.set_yscale('log')
    # save
    fig.savefig(home+'power_{}.png'.format(wlim),bbox_inches='tight')
