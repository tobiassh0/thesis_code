

def plot_example(home,simloc,wkmax=[10,20],minspec='Protons'):
    # load FT2d of one case (5%)
    sim = getSimulation(home+simloc)
    d0 = sdfread(0)
    times = read_pkl('times')
    FT2d = read_pkl('FT_2d_Magnetic_Field_Bz')
    Nw, Nk = FT2d.shape
    klim, wlim = getDispersionlimits(sim)
    kmin, kmax = klim
    wmin, wmax = wlim
    # frequency and Alfven speed
    wcmin = getCyclotronFreq(d0,minspec)
    vA = getAlfvenVel(d0)
    # normalised cell widths
    dk = (2*const.PI/getGridlen(d0))
    dw = (2*const.PI/times[-1])
    # extents = [0,kmax*vA/wcmin,0,wmax/wcmin]
    # thresh FT2d
    wmax_prime, kmax_prime = wkmax
    extents = [0,kmax_prime,0,wmax_prime]
    # threshold FT2d
    kthresh = int(FT2d.shape[1]*extents[1]*wcmin/vA/kmax)
    wthresh = int(FT2d.shape[0]*extents[-1]*wcmin/wmax)
    tFT2d = FT2d[:wthresh,:kthresh]
    tNw, tNk = tFT2d.shape
    kx = np.linspace(0,kmax_prime,100)

    # setup gridspec
    fig = plt.figure(figsize=(12,5))
    # FT2d imshow
    gs1 = GridSpec(4, 6, bottom=0.02, top=0.98, left=0.02, right=0.47, wspace=0.0, hspace=0.0)# bottom spacing creates space for gs2 
    ax1 = fig.add_subplot(gs1[:, :])
    ax1.imshow(np.log10(tFT2d[:,1:]),**imkwargs,extent=extents,cmap='magma')
    # power along line
    gs2 = GridSpec(8, 6, bottom=0.02, top=0.98, left=0.55, right=0.98, wspace=0.25, hspace=0.25) # nrows, ncols, l, r, wsp, hsp
    ax4 = fig.add_subplot(gs2[4:, :3])
    ax5 = fig.add_subplot(gs2[4:, 3:])
    ax2 = fig.add_subplot(gs2[:4, :3])
    ax3 = fig.add_subplot(gs2[:4, 3:])
    axarr = [ax2,ax3,ax4,ax5]
    labels= ['(a)','(b)','(c)','(d)']

    # line doppler on 4 example lines (same color and style)
    freqs = np.linspace(0,wmax_prime,tNw)
    ynarr = tNw*freqs/wmax_prime
    angles = [-np.pi/2,-np.pi/2.5,-np.pi/3,-np.pi/4] # -90, -72, -60, -45
    power_lines = np.zeros((len(angles),len(freqs)))
    # plot all cases
    for i in range(len(angles)):
        for j in range(0,wmax_prime+1):
            axarr[i].axvline(j,linestyle='--',color='darkgrey')
        power_lines[i,:] = pl.getPowerLine(tFT2d,ynarr,freqs,kmax_prime=100,wmax_prime=wmax_prime,angle=angles[i])
        axarr[i].plot(freqs,np.log10(power_lines[i,:])-np.log10(dw),color='k')
        axarr[i].set_ylim(-3,0.5) ; axarr[i].set_xlim(0,wmax_prime)
        ax1.plot(kx,kx*(1/np.tan(angles[i]))+8,color='k',linestyle='--')
        if i==0: 
            y0=8.125 ; x0=18
        elif i==1:
            y0=2.225 ; x0=18
        else:
            y0=0.125 ; x0=(-8*np.tan(angles[i]))
        ax1.annotate(labels[i],xy=(x0,y0),xycoords='data',**tnrfont)
        axarr[i].annotate(labels[i],xy=(0.05,0.95),xycoords='axes fraction',va='top',ha='left',**tnrfont)
    # # plot actual power (using kmax_prime~=1000)
    # log10_power = read_pkl('log10_power')
    # omegas = read_pkl('omegas_power')/wcmin
    # plt.plot(omegas,log10_power-np.log10(dw),color='r')

    # formatting 
    ax1.set_xlim(0,kmax_prime)
    ax1.set_ylim(0,wmax_prime)
    ax2.set_xticklabels([])
    ax3.set_xticklabels([])
    ax3.set_yticklabels([])
    ax5.set_yticklabels([])
    ax1.set_xlabel(r'$kv_A/\Omega_p$',**tnrfont)
    ax1.set_ylabel(r'$\omega/\Omega_p$',**tnrfont)
    ax4.set_ylabel('PSD',y=1,**tnrfont)
    ax4.set_xlabel(r'$\omega^\prime/\Omega_p$',x=1,**tnrfont)
    plt.savefig(home+'FT2d_line_example.png',bbox_inches='tight')
    # plt.show()
    pass


if __name__=='__main__':
    from func_load import *
    import dopplerShift.line_integrate as li
    import power.power_line as pl
    from matplotlib.gridspec import GridSpec


    plot_example(home='/storage/space2/phrmsf/lowres_D_He3/',simloc='0_05_p_90')
