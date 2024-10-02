
from func_load import *

if __name__=='__main__':
    home = '/storage/space2/phrmsf/lowres_D_He3/'
    sims = [i for i in os.listdir(home) if 'p_90' in i]
    mD = getMass('Deuterons') ; mHe3 = getMass('He3')
    qqD = (const.qe)**2 ; qqHe3 = (2*const.qe)**2
    for sim in sims:
        print(sim)
        simloc=getSimulation(home+sim)
        d0=sdfread(0)
        # load times
        times=read_pkl('times')
        tcp=2*const.PI/getCyclotronFreq(d0,'Protons')
        # Alfven velocity
        vA=getAlfvenVel(d0)
        # get all species
        all_species=getIonSpecies(d0)[2:] # ignore minority species
        # find which files are restart files
        rest_files=para_check_restart(simloc)
        rest_times=[]
        # initialise empty species arrays
        Vxpara = np.zeros((2,len(rest_files)))
        Vypara = np.zeros((2,len(rest_files)))
        Vzpara = np.zeros((2,len(rest_files)))
        # get magnetic field angle
        theta_xz, theta_xy = getMagneticAngle(d0)
        # loop through all restart files and species
        for i in range(len(rest_files)):
            d=sdfread(rest_files[i])
            rest_times.append(d.__dict__['Header']['time'])
            for k in range(len(all_species)):
                _, [vxpara, vypara, vzpara] = getPerpParaVel(d,all_species[k],theta_xz=theta_xz,theta_xy=theta_xy)
                # print(vxpara, vypara, vzpara)
                Vxpara[k,i] = vxpara
                Vypara[k,i] = vypara
                Vzpara[k,i] = vzpara
        # plot vpara through time
        fig,ax=plt.subplots(figsize=(8,6))
        # ax.plot(rest_times,Vxpara[0,:],'bo-')
        # ax.plot(rest_times,Vypara[0,:],'ro-')
        # ax.plot(rest_times,Vzpara[0,:],'go-')
        ax.plot(rest_times/tcp,np.sqrt(Vxpara[0,:]**2 + Vypara[0,:]**2 + Vzpara[0,:]**2)/vA,'ko-')
        ax.set_ylabel(r'$|v_\parallel|/v_A$',**tnrfont)
        ax.set_xlabel('Time, '+r'$\tau_{cp}$',**tnrfont)
        fig.savefig(home+sim+'/vpara_time.png')
        plt.clf()