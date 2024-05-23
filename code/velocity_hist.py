
from func_load import *
import energy.LarmorRadii as lr

def ciggie_vperp(restart_files,vA,spec,bins=100,theta=3.141/2,binrange=(0,0.2)):
    # create empty array of size len(restart times) x No. bins
    cig_mat = np.zeros((len(restart_files),bins))
    # get mass
    mass_spec = getMass(spec)
    # loop through times
    for i in range(len(restart_files)):
        print(i,restart_files[i])
        di = sdfread(restart_files[i])
        # load velocities
        Vx = getQuantity1d(di,'Particles_Px_'+spec)/mass_spec # getQuantity1d
        Vy = getQuantity1d(di,'Particles_Py_'+spec)/mass_spec
        Vz = getQuantity1d(di,'Particles_Pz_'+spec)/mass_spec
        # perp
        Vxperp = Vx*np.sin(theta)           # (1-(np.cos(theta)**2)*(np.cos(theta_y)**2))**0.5
        Vyperp = Vy                         # (1-(np.cos(theta)**2)*(np.sin(theta_y)**2))**0.5
        Vzperp = Vz*np.cos(theta)
        Vperp = np.sqrt(np.mean(Vxperp**2 + Vyperp**2 + Vzperp**2))
        print(Vperp)
        cig_mat[i,:],_ = np.histogram(Vperp/vA,bins=bins,density=True,range=binrange)
    return cig_mat

if __name__=='__main__':

    simloc = getSimulation('/storage/space2/phrmsf/lowres_D_He3/0_05_p_90')
    times = read_pkl('times')
    restart_files = lr.para_check_restart(simloc)
    # restart_files = [0]#,3000,6000,9000,12000]
    d0 = sdfread(0)
    B0 = getMeanField3D(d0,'Magnetic_Field_B')
    theta,_ = getMagneticAngle(d0) # radians
    vA = getAlfvenVel(d0)

    binrange=(0,0.1)
    cig_mat = ciggie_vperp(restart_files,vA,spec='Deuterons',bins=1000,theta=theta,binrange=binrange)
    plt.imshow((cig_mat.T),**imkwargs,cmap='jet',extent=[0,10,binrange[0],binrange[1]])
    plt.show()
    sys.exit()

    # loop through and get histogram of y-velocity
    # N,M=(5,5)
    # fig,ax=plt.subplots(figsize=(15,15),nrows=N,ncols=M,sharex=True,sharey=True)
    # ax=ax.ravel()
    for i in range(len(restart_files)):
        print(i)
        rstfl = sdfread(restart_files[i])
        DVx  =np.sin(theta)*getQuantity1d(rstfl,'Particles_Px_Deuterons')/getMass('Deuterons')
        DVy  =getQuantity1d(rstfl,'Particles_Py_Deuterons')/getMass('Deuterons')
        DVz  =np.cos(theta)*getQuantity1d(rstfl,'Particles_Pz_Deuterons')/getMass('Deuterons')
        He3Vx=np.sin(theta)*getQuantity1d(rstfl,'Particles_Px_He3')/getMass('He3')
        He3Vy=getQuantity1d(rstfl,'Particles_Py_He3')/getMass('He3')
        He3Vz=np.cos(theta)*getQuantity1d(rstfl,'Particles_Pz_He3')/getMass('He3')
        _=plt.hist(DVx/vA,bins=100,color='b',alpha=0.45,density=True,range=(-0.2,0.2))
        _=plt.hist(He3Vx/vA,bins=100,color='r',alpha=0.45,density=True,range=(-0.2,0.2))
        _=plt.hist(He3Vx**2)
    # fig.supylabel('Normalised Count',**tnrfont)
    # fig.supxlabel(r'$v_y/v_A$',**tnrfont)
    plt.show()

