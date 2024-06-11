from func_load import *
import multiprocessing as mp # for parallelisation
import matplotlib.cm as cm

def para_loop_restart(i,check_species='Protons'):
    if 'Particles_Px_'+check_species in getKeys(sdfread(i)):
        return i
    else:
        pass

def para_check_restart(simloc):
    print('Parallel checking restart files...')
    pool=mp.Pool(mp.cpu_count())
    rest_files = np.array(pool.map_async(para_loop_restart,list_sdf(simloc)).get(99999))
    pool.close()
    # remove None
    rest_files = rest_files[rest_files != np.array(None)]
    print('Done :: {} files found.\n{}'.format(len(rest_files),rest_files))
    return rest_files

def linear_check_getVelocityPerp(restart_files,times,species,theta,theta_y=0):
    print('Linear loading velocities of {} for {} files...'.format(species,len(restart_files)))
    Vperp = np.zeros(len(restart_files))
    Vperp_squared = []
    restart_times = np.zeros(len(restart_files))
    mass_spec = getMass(species)
    for i in range(len(restart_files)):
        restart_times[i] = times[restart_files[i]]
        di = sdfread(restart_files[i])
        print(i,restart_files[i])
        # load velocities
        Vx = getQuantity1d(di,'Particles_Px_'+species)/mass_spec # getQuantity1d
        Vy = getQuantity1d(di,'Particles_Py_'+species)/mass_spec
        Vz = getQuantity1d(di,'Particles_Pz_'+species)/mass_spec
        # perp
        Vxperp = Vx*np.sin(theta)           # (1-(np.cos(theta)**2)*(np.cos(theta_y)**2))**0.5
        Vyperp = Vy                         # (1-(np.cos(theta)**2)*(np.sin(theta_y)**2))**0.5
        Vzperp = Vz*np.cos(theta)           # this is unchanged following angle in xy plane (cylindrical coords)
        # # para
        # Vxpara = Vx*np.cos(theta)           # (1-(np.cos(theta)**2)*(np.cos(theta_y)**2))**0.5
        # Vypara = 0                         # (1-(np.cos(theta)**2)*(np.sin(theta_y)**2))**0.5
        # Vzpara = Vz*np.sin(theta)
        Vperp_squared.append(Vxperp**2 + Vyperp**2 + Vzperp**2)
        # Vperp[i] = np.sqrt(np.mean(Vxperp**2 + Vyperp**2 + Vzperp**2))

    return np.array(Vperp_squared), restart_times

def plotChangeEnergyLarmorRatio(home,sims,species,minspec='Protons',labels=[''],colors=['k']):
    if len(colors) != len(sims):
        colors=cm.rainbow(np.linspace(0,1,len(sims)))

    for i in range(len(sims)):
        simloc = getSimulation(home+sims[i])
        times = read_pkl('times')
        tcmin = 2*const.PI/getCyclotronFreq(sdfread(0),minspec)
        masses = [getMass(i) for i in species]
        charges = [getChargeNum(i)*const.qe for i in species]
        uarr = [read_pkl(i+'_KEdens') for i in species]
        xiarr = getConcentrationRatios(sdfread(0))
        Earr = [uarr[i]/xiarr[i] for i in range(len(species))]
        dEarr = [(Earr[i]-np.mean(Earr[i][::10])) for i in range(len(species))]
        drratio = (dEarr[0]/dEarr[1])*(masses[0]/masses[1])*(charges[1]/charges[0])**2
        plt.plot(times/tcmin,drratio,color=colors[i],label=labels[i])
    plt.yscale('symlog')
    plt.xlabel(r'$t/\tau_{cp}$',**tnrfont)
    plt.ylabel(r'$(\Delta E_D/\Delta E_{He3})(m_D/m_{He3})(q_{He3}/q_D)^2$',**tnrfont)
    plt.legend(loc='best')
    plt.xlim(0,times[-1]/tcmin)
    plt.savefig(home+'EqualChangeLarmorRadii.png',bbox_inches='tight')
    plt.show()
    return None

def plotChangeLarmorRadius(home,sims,species=['Deuterons','He3'],minspec='Protons',labels=[''],colors=['b','r']):
    for i in range(len(sims)):
        simloc=getSimulation(home+sims[i])
        d0=sdfread(0)
        times=read_pkl('times')
        tcmin=2*const.PI/getCyclotronFreq(d0,minspec)
        # find restart files where Vx, Vy and Vz are present # loop through all sdf files, see if Vx_Particle in keys
        restart_files = para_check_restart(simloc)

        # get magnetic field angle (wrt y) and field strength # TODO; generalise for B field rotated wrt y and x
        theta,theta_y=getMagneticAngle(d0)
        B0 = getMeanField3D(d0,'Magnetic_Field_B')
        LDe = getDebyeLength(d0, 'Electrons')

        # linear loading velocities
        for j in range(len(species)):
            mass_spec = getMass(species[j])
            # get Vperp and times at restart
            Vperp_squared,restart_times = linear_check_getVelocityPerp(restart_files,times,species[j],theta,theta_y)
            
            # # 2d Larmor array
            # rLarr = np.sqrt(Vperp_squared)*mass_spec/(getChargeNum(species[j])*const.qe*B0)
            # plt.imshow(rLarr,**imkwargs)
            
            # RMS Vperp
            Vperp = np.sqrt(np.mean(Vperp_squared,axis=1)) # mean Vperp per time step
            
            # plot rL^2 as function of time & total energy transferred (should be equal)
            rL2 = (mass_spec*Vperp/(getChargeNum(species[j])*const.qe*B0))**2
            print(rL2.shape,restart_times.shape)
            plt.plot(restart_times/tcmin,(rL2-rL2[0])/LDe**2,color=colors[j],marker='o')

        # format and save
        plt.xlabel(r'$t/\tau_{cp}$',**tnrfont)
        plt.ylabel(r'$\Delta r^2_{L\sigma}/\lambda_{De}$',**tnrfont)
        plt.title(labels[i]+"%",**tnrfont)
        plt.legend(species,loc='best')
        plt.savefig(home+'Delta_LarmorRadiiSquared_'+labels[i]+'.png')
    return None

if __name__=='__main__':

    home = '/storage/space2/phrmsf/lowres_D_He3/'
    # load sim
    xiHe3arr = ['05']
    xiHe3arr = ['05','10','15','22','25','34','38','45']
    sims = ['0_'+xiHe3+'_p_90' for xiHe3 in xiHe3arr]
    species = ['Deuterons','He3']
    
    # # change in kinetic energy assuming equal change in larmor radii
    # plotChangeEnergyLarmorRatio(home,sims,species,labels=xiHe3arr)
    # sys.exit()

    # plot change in Larmor radius
    plotChangeLarmorRadius(home,sims,species,labels=xiHe3arr)
    sys.exit()
