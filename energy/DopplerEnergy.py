
def getDopplerEnergy(sims,Eminarr,wkmax=[10,20],minspec='Protons',quantity='Magnetic_Field_Bz',vperp_vA=0.9):
    fig,axs=plt.subplots(figsize=(10,10),nrows=2,ncols=2,sharex=True,sharey=True,layout='constrained')
    axs=axs.ravel()
    for i in range(len(sims)):
        # load sim
        _ = getSimulation(sims[i])
        times = read_pkl('times')
        d0 = sdfread(0)
        dk = (2*const.PI/getGridlen(d0))
        dw = (2*const.PI/times[-1])
        vA = getAlfvenVel(d0)
        wcmin = getCyclotronFreq(d0,minspec)
        klim = 0.5*2*const.PI/getdxyz(d0)
        wlim = 0.5*2*const.PI/getdt(times)
        FT2d = read_pkl('FT_2d_'+quantity)
        # thresh FT2d
        (nw,nk) = FT2d.shape
        wmax = wkmax[0]*wcmin ; kmax = wkmax[1]*wcmin/vA
        tFT2d = FT2d[:int(nw*wmax/wlim),:int(nk*kmax/klim)]
        # tFT2d[np.log10(tFT2d) < 1.2] = 0
        (tnw,tnk) = tFT2d.shape
        im=axs[i].imshow(np.log10(tFT2d[:,1:]),**imkwargs,cmap='magma',extent=[0,wkmax[1],0,wkmax[0]],clim=(1.2,4))
        # # line integrate gradient
        # grad, angle = li.getdopgrad(tFT2d,logthresh=1.8,norm=(dk*vA/dw),kernel='custom',dydxrange=(-0.5,0.5)) 

        # calc theory Dopp vel
        umin = np.sqrt((1e6*Eminarr[i]*const.qe)/(0.5*getMass(minspec)))
        print(vA,wcmin)

        # pitch-angle
        theta_B,_ = getMagneticAngle(d0)
        phi = np.arcsin((vperp_vA*vA)/umin)
        vdop = (umin/vA)*np.cos(theta_B)*np.cos(phi)
        print(vdop)

        kx = np.linspace(0,wkmax[1],100)*wcmin/vA
        # plot threshed FT2d and Dopp lines
        for n in range(0,wkmax[0]+3):
            w_prime = n*wcmin - kx*np.abs(vdop*vA)
            axs[i].plot(kx*vA/wcmin,w_prime/wcmin,color='w',linestyle='--')
        # set title
        axs[i].set_title(str(Eminarr[i])+'MeV',**tnrfont)

    axs[0].set_xlim(0,wkmax[1]) ; axs[0].set_ylim(0,wkmax[0])
    fig.supylabel(r'$\omega/\Omega_p$',**tnrfont)
    fig.supxlabel(r'$kv_A/\Omega_p$',**tnrfont)

    # colorbar
    p0 = axs[0].get_position().get_points().flatten()
    p1 = axs[1].get_position().get_points().flatten()
    p3 = axs[-1].get_position().get_points().flatten()
    ax_cbar = fig.add_axes([1, 0.0725, 0.02, 0.895]) # [left bottom width height]
    cbar = plt.colorbar(im, cax=ax_cbar, orientation='vertical')
    # plt.show()
    return fig

if __name__=='__main__':
    from func_load import *
    import dopplerShift.getDopplerShiftICE as gDs
    import dopplerShift.line_integrate as li
    import os,sys

    home = '/storage/space2/phrmsf/lowres_D_He3/energy_protons/'
    sims = [i for i in os.listdir(home) if 'MEV' in i]
    proton_energies = [int(i.split('M')[0]) for i in sims]
    # sort in order of increasing energy
    proton_energies, sims = zip(*sorted(zip(proton_energies,sims)))
    print(sims,proton_energies)

    # plot Doppler
    fig = getDopplerEnergy([home+s for s in sims][1:],proton_energies[1:])
    fig.savefig(home+'FT2d_Doppler.png',bbox_inches='tight')
