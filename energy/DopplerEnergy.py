
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
    print(sims)

    # getdoptheory(file0,Emin,wcmin,minspec,vA,uperp_vA=0.9)

    FT2darr = []
    grads = []
    quantity = 'Magnetic_Field_Bz'
    for i in range(len(sims)):
        # get Doppler velocity
        _ = getSimulation(home+sims[i])
        times = read_pkl('times')
        d0 = sdfread(0)
        dk = (2*const.PI/getGridlen(d0))
        dw = (2*const.PI/times[-1])
        vA = getAlfvenVel(d0)
        wcmin = getCyclotronFreq(d0,'Protons')
        klim = 0.5*2*const.PI/getdxyz(d0)
        wlim = 0.5*2*const.PI/getdt(times)
        FT2d = read_pkl('FT_2d_'+quantity)
        FT2darr.append(FT2d)
        # thresh FT2d
        (nw,nk) = FT2d.shape
        wmax = 10*wcmin ; kmax = 20*wcmin/vA
        tFT2d = np.log10(FT2d[:int(nw*wmax/wlim),:int(nk*kmax/klim)])
        (tnw,tnk) = tFT2d.shape
        grad, angle = li.getdopgrad(tFT2d,logthresh=1.8,norm=(dk*vA/dw),kernel='custom',dydxrange=(-0.5,0.5)) 
        grads.append(grad)
    # sims = [home+i for i in sims]
    # _=gDs.getKernelDoppler(sims,FT2darr,normspecies='Protons',plot=False,home=home,logthresh=1.5)
    # dsv_vA,dsv,vA = dsvarr[:,0], dsvarr[:,1], dsvarr[:,2]

    # plot as a function of proton energy
    print(grads)
    plt.scatter(proton_energies,grads)
    plt.xlabel('Energy [MeV]',**tnrfont)
    plt.ylabel('vdop [vA]',**tnrfont)
    plt.show()
