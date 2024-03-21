
if __name__=='__main__':
    from func_load import *
    from dopplerShift.getDopplerShiftICE import *
    import os,sys

    home = '/storage/space2/phrmsf/lowres_D_He3/energy_protons/'
    sims = [i for i in os.listdir(home) if 'MEV' in i]
    proton_energies = [int(i.split('M')[0]) for i in sims]
    # sort in order of increasing energy
    proton_energies, sims = zip(*sorted(zip(proton_energies,sims)))
    print(sims)
    FT2darr = [] ; quantity = 'Magnetic_Field_Bz'
    for i in range(len(sims)):
        # get Doppler velocity
        getSimulation(home+sims[i])
        FT2darr.append(read_pkl('FT_2d_'+quantity))
    dsvarr = getKernelDoppler(sims,FT2darr,normspecies='Protons',plot=False,home=home,logthresh=1.5)
    # plot as a function of proton energy
    print(dsvarr)
    dsv_vA,dsv,vA = dsvarr[:,0], dsvarr[:,1], dsvarr[:,2]
    plt.plot(proton_energies,dsv)
    plt.xlabel('Energy [MeV]',**tnrfont)
    plt.ylabel('dsv [m/s]',**tnrfont)
    plt.show()
