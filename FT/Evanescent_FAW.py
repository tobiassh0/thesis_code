
from func_load import *

def getEvanescentMCI(home,sim,wklim=[30,100],species=['Electrons','Deuterons','Alphas'],theta=None,plot=False):
    figFT2d, axFT2d = plt.subplots(figsize=(8,6))
    if theta==None:
        theta = 89*const.PI/180
    print(sim)
    simloc = getSimulation(home+sim)
    FT2d = read_pkl('FT_2d_Magnetic_Field_Bz')
    d0=sdfread(0)
    times=read_pkl('times')
    knyq = 0.5*2*const.PI/getdxyz(d0)
    wnyq = 0.5*2*const.PI/getdt(times)
    vA = getAlfvenVel(d0)
    wcD = getCyclotronFreq(d0,'Deuterons')
    tFT2d = FT2d[:int(FT2d.shape[0]*(wklim[0]*wcD)/wnyq),:int(FT2d.shape[1]*(wklim[1]*wcD/vA)/knyq)]
    del FT2d
    axFT2d.imshow(np.log10(tFT2d),**imkwargs,extent=[0,wklim[1],0,wklim[0]],cmap='magma')
    
    # cold plasma dispersion (complex k)
    omegas=np.linspace(0,wklim[0]*wcD,10000)
    wpf=[getPlasmaFreq(d0,i) for i in species]
    wcf=[getCyclotronFreq(d0,i) for i in species]
    print(wpf,wcf)
    k1,k2,k3=coldplasmadispersion_analytical(omegas,wpf=wpf,wcf=wcf,theta=theta)
    # # plt.scatter(np.real(k1)*vA/wcD,omegas/wcD,color='b')
    # # plt.scatter(np.imag(k2)*vA/wcD,omegas/wcD,color='r')
    # # plt.scatter(np.real(k3)*vA/wcD,omegas/wcD,color='b')
    # # plt.scatter(np.imag(k3)*vA/wcD,omegas/wcD,color='r')
    if plot:
        axFT2d.plot(np.imag(k1)*vA/wcD,omegas/wcD,color='r')
        axFT2d.plot(np.real(k2)*vA/wcD,omegas/wcD,color='b')

    # find positions (above some thresh freq.) of phase matching
    thresh = omegas/wcD > 10
    index = np.argmin(np.abs(np.imag(k1[thresh]) - np.real(k2[thresh])))
    print(omegas[thresh][index]/wcD)
    if plot:
        axFT2d.set_xlabel(r'$kv_A/\Omega_D$',**tnrfont)
        axFT2d.set_ylabel(r'$\omega/Omega_D$',**tnrfont)
        axFT2d.scatter(k2[thresh][index]*vA/wcD,omegas[thresh][index]/wcD,color='k',zorder=0)
        axFT2d.set_xlim(0,wklim[1])
        axFT2d.set_ylim(0,wklim[0])
        figFT2d.savefig(home+sim+'/FT2d_MCI_Evanescent.png',bbox_inches='tight')
    return omegas[thresh][index]/wcD

def plotEvanescentMCI(home,sims,xi2):
    omega_Ev=np.zeros(len(sims))
    for i in range(len(sims)):
        omega_Ev[i]=getEvanescentMCI(home,sims[i],plot=True)
    plt.clf()
    plt.plot(xi2,omega_Ev,'ro-')
    plt.xlabel(r'$\xi_2$',**tnrfont)
    plt.ylabel(r'$\omega_{ev.}/\Omega_D$',**tnrfont)
    plt.savefig(home+'Evansecent_freqs.png',bbox_inches='tight')
    return None

def getPowerEvanescentMCI(home,sim):
    simloc = getSimulation(home+sim)
    d0=sdfread(0)
    wcD=getCyclotronFreq(d0,'Deuterons')
    times=read_pkl('times')
    wnyq=0.5*2*const.PI/getdt(times)
    dw=2*wnyq/len(times)
    figPSD, axPSD = plt.subplots(figsize=(8,4))
    # vertical harmonic lines
    for i in range(0,31):
        axPSD.axvline(i,color='darkgrey',linestyle='--')
    # power
    omegas_power=read_pkl('omegas_power')
    log10_power=read_pkl('log10_power')
    axPSD.plot(omegas_power/wcD,np.log10((10**log10_power)/dw),color='k')
    # Evanescent freq.
    omega_Ev=getEvanescentMCI(home,sim)
    axPSD.axvline(omega_Ev,color='red',linestyle='--')
    axPSD.set_xlim(0,30)
    figPSD.savefig(home+sim+'/Evanescent_freq_PSD.png',bbox_inches='tight')
    return None

def plotPowerEvanescentMCI(home,sims):
    plt.clf()
    for i in range(len(sims)):
        getPowerEvanescentMCI(home,sims[i])
    return None

if __name__=='__main__':
    home = '/storage/space2/phrmsf/traceT/'
    disallowed = ['traceT_D_100_T_00_v3']
    sims = [i for i in os.listdir(home) if 'traceT_D_' in i and i not in disallowed]
    xiT = [int(i.split('_')[4]) for i in sims]
    xiT, sims = sort_arrays(xiT,sims)
    plotEvanescentMCI(home,sims,xiT)
    plotPowerEvanescentMCI(home,sims)

    sys.exit()
    for sim in sims:
        simloc = getSimulation(home+sim)
        print(sim)
        FT2d = read_pkl('FT_2d_Magnetic_Field_Bz')
        d0=sdfread(0)
        times=read_pkl('times')
        knyq = 0.5*2*const.PI/getdxyz(d0)
        wnyq = 0.5*2*const.PI/getdt(times)
        dw = wnyq/(2*len(times))
        vA = getAlfvenVel(d0)
        wcD = getCyclotronFreq(d0,'Deuterons')
        wklim = [30,100] # normalised plotting limits
        tFT2d = FT2d[:int(FT2d.shape[0]*(wklim[0]*wcD)/wnyq),:int(FT2d.shape[1]*(wklim[1]*wcD/vA)/knyq)]
        del FT2d
        plt.imshow(np.log10(tFT2d),**imkwargs,extent=[0,wklim[1],0,wklim[0]],cmap='magma')
        
        # cold plasma dispersion (complex k)
        theta=89*const.PI/180
        species=['Electrons','Deuterons','Alphas']
        omegas=np.linspace(0,wklim[0]*wcD,10000)
        wpf=[getPlasmaFreq(d0,i) for i in species]
        wcf=[getCyclotronFreq(d0,i) for i in species]
        print(wpf,wcf)
        k1,k2,k3=coldplasmadispersion_analytical(omegas,wpf=wpf,wcf=wcf,theta=theta)
        # # plt.scatter(np.real(k1)*vA/wcD,omegas/wcD,color='b')
        # plt.scatter(np.imag(k1)*vA/wcD,omegas/wcD,color='r')
        # plt.scatter(np.real(k2)*vA/wcD,omegas/wcD,color='b')
        # # plt.scatter(np.imag(k2)*vA/wcD,omegas/wcD,color='r')
        # # plt.scatter(np.real(k3)*vA/wcD,omegas/wcD,color='b')
        # # plt.scatter(np.imag(k3)*vA/wcD,omegas/wcD,color='r')

        # find positions (above some thresh freq.) of phase matching
        thresh = omegas/wcD > 10
        index = np.argmin(np.abs(np.imag(k1[thresh]) - np.real(k2[thresh])))
        # print(omegas[thresh][index]/wcD)
        # plt.scatter(k2[thresh][index]*vA/wcD,omegas[thresh][index]/wcD,color='k')
        # plt.ylim(0,wklim[0])
        # plt.xlim(0,wklim[1])
        # plt.show()
        plt.scatter(xiT[i],omegas[thresh][index]/wcD,color='r')
        i+=1
    plt.show()
