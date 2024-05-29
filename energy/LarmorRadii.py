from func_load import *
import multiprocessing as mp # for parallelisation

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
    print('Done :: {} files found.'.format(len(rest_files)))
    return rest_files

def linear_check_getVelocityPerp(restart_files,times,spec,mass_spec,theta,theta_y=3.1415/2):
    print('Linear loading velocities of {} for {} files...'.format(spec,len(restart_files)))
    Vperp = np.zeros(len(restart_files))
    Vperp_squared = []
    restart_times = np.zeros(len(restart_files))
    for i in range(len(restart_files)):
        restart_times[i] = times[restart_files[i]]
        di = sdfread(restart_files[i])
        print(i,restart_files[i])
        # load velocities
        Vx = getQuantity1d(di,'Particles_Px_'+spec)/mass_spec # getQuantity1d
        Vy = getQuantity1d(di,'Particles_Py_'+spec)/mass_spec
        Vz = getQuantity1d(di,'Particles_Pz_'+spec)/mass_spec
        # perp
        Vxperp = Vx*np.sin(theta)           # (1-(np.cos(theta)**2)*(np.cos(theta_y)**2))**0.5
        Vyperp = Vy                         # (1-(np.cos(theta)**2)*(np.sin(theta_y)**2))**0.5
        Vzperp = Vz*np.cos(theta)
        # # para
        # Vxpara = Vx*np.cos(theta)           # (1-(np.cos(theta)**2)*(np.cos(theta_y)**2))**0.5
        # Vypara = 0                         # (1-(np.cos(theta)**2)*(np.sin(theta_y)**2))**0.5
        # Vzpara = Vz*np.sin(theta)
        Vperp_squared.append(Vxperp**2 + Vyperp**2 + Vzperp**2)
        # Vperp[i] = np.sqrt(np.mean(Vxperp**2 + Vyperp**2 + Vzperp**2))

    return np.array(Vperp_squared), restart_times


if __name__=='__main__':
    home = '/storage/space2/phrmsf/lowres_D_He3/'
    # load sim
    xiHe3arr = ['05','10','15','22','25','34','38','45']
    xiHe3arr = ['05']
    for xiHe3 in xiHe3arr:
        simloc=getSimulation(home+'0_'+xiHe3+'_p_90')
        d0=sdfread(0)
        times=read_pkl('times')
        tcmin=2*const.PI/getCyclotronFreq(d0,'Protons')
        # find restart files where Vx, Vy and Vz are present # loop through all sdf files, see if Vx_Particle in keys
        # for fl in list_sdf(simloc):
        #     if 'Particles_Px_'+spec in getKeys(sdfread(fl)):
        #         rest_files.append(fl)
        restart_files = para_check_restart(simloc)

        # get magnetic field angle (wrt y) and field strength # TODO; generalise for B field rotated wrt y and x
        theta,theta_y=getMagneticAngle(d0)
        B0 = getMeanField3D(d0,'Magnetic_Field_B')
        LDe = getDebyeLength(d0, 'Electrons')

        # get Velocities and Larmor Radii
        species = ['Deuterons','He3']
        colors = ['b','r']
        j=0
        # linear loading velocities
        for spec in species:
            mass_spec = getMass(spec)
            # get Vperp and times at restart
            Vperp_squared,restart_times = linear_check_getVelocityPerp(restart_files,times,spec,mass_spec,theta,theta_y)
            # 2d Larmor array
            rLarr = np.sqrt(Vperp_squared)*mass_spec/(getChargeNum(spec)*const.qe*B0)
            plt.imshow(rLarr,**imkwargs)
            # # RMS Vperp
            # Vperp = np.sqrt(np.mean(Vperp_squared,axis=1)) # mean Vperp per time step
            # # plot rL as function of time & total energy transferred
            # rL = mass_spec*Vperp/(getChargeNum(spec)*const.qe*B0)
            # print(rL.shape,restart_times.shape)
            # plt.plot(restart_times/tcmin,(rL-rL[0])/LDe,color=colors[j],marker='o')
            j+=1
        plt.ylabel('Normalised Count',**tnrfont)
        plt.xlabel(r'$r_{L\sigma}/\lambda_{De}$',**tnrfont)
        plt.show()
        # # show
        # plt.xlabel(r'$t/\tau_{cp}$',**tnrfont)
        # plt.ylabel(r'$\Delta r_{L\sigma}/\lambda_{De}$',**tnrfont)
        # plt.title(xiHe3+"%",**tnrfont)
        # plt.legend(species,loc='best')
        # plt.savefig(home+'Delta_LarmorRadii_'+xiHe3+'.png')
        # plt.clf()
