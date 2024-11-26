
from func_load import * 

class tS_Simulation(object):
    def __init__(self,home,loc):
        self.home = home
        self.simloc = loc

        # all files
        self.ind_lst = list_sdf(self.simloc)
        self.d0 = sdfread(0)
        self.times = read_pkl('times')
        #self.tS_gettimes()
        dt = (self.times[-1]/len(self.times))
        
        # for k in getKeys(d0): print(k)
        # raise SystemExit

        # Nyquist, velocities and normalisations
        self.wpe = getPlasmaFreq(self.d0,'Left_Electrons')
        self.tpe=2*const.PI/self.wpe
        self.u0x=np.abs(getMeanquantity(self.d0,'Particles_Vx_Left_Electrons'))
        self.kmax=0.5*2*const.PI/getdxyz(self.d0) # nyquist wavenum.
        self.kmax_prime=self.kmax*self.u0x/self.wpe
        self.wmax=0.5*2*const.PI/dt # nyquist freq.
        self.wmax_prime=self.wmax/self.wpe

        # axes labels
        self.tlabel = r'$t\omega_{pe}/2\pi$'
        self.klabel = r'$ku_{0x}/\omega_{pe}$'
        self.wlabel = r'$\omega/\omega_{pe}$'
    
    # get times
    def tS_gettimes(self):
        times=[]
        for i in range(len(self.ind_lst)):
            times.append(sdfread(i).__dict__['Header']['time'])
        return times

    # load fields
    def tS_LoadEX(self):
        _,tEx = get_batch_fieldmatrix(self.ind_lst,['Electric_Field_Ex'],'Electric_Field_Ex',load=True,para=False)
        return tEx

    # load Bz field
    def tS_LoadBZ(self):
        _,tBz = get_batch_fieldmatrix(self.ind_lst,load=True,para=False)
        return tBz

    # 1d FT
    def FT1d(self):
        # load Ex ::
        Ex = self.tS_LoadEX()
        # make ::
        FT1d = get1dTransform(Ex,window=False)
        # plot ::
        fig1d,ax1d=plt.subplots(figsize=(6,6))
        im=ax1d.imshow(np.log10(FT1d[1:,:]),**imkwargs,cmap='magma',\
                        extent=[0,self.kmax_prime,0,self.times[-1]/self.tpe])
        # cbar=plt.colorbar(im)
        ax1d.set_xlabel(self.klabel,**tnrfont)
        ax1d.set_ylabel(self.tlabel,**tnrfont)
        ax1d.set_xlim(0,15)
        ax1d.axvline(3**0.5/2,linestyle='--',color='k') # kmax from theory
        fig1d.savefig(self.home+'/twoStream_FT1d_Ex.png',bbox_inches='tight')
        plt.show()
    # 2d FT
    def FT2d(self):
        # load Bz ::
        Bz = self.tS_LoadBZ()
        # make ::
        tFT2d = np.fft.fft2(HanningWindowT(Bz))[:,:] # windowed in time # windowed in both # HanningWindow2D
        tFT2d = np.fft.fftshift(tFT2d)
        tFT2d = tFT2d[int(Bz.shape[0]/2):,:] # k > 0
        FT2d = np.abs(tFT2d) # destroy phase information
        del tFT2d
        # plot ::
        fig2d,ax2d=plt.subplots(figsize=(6,6),layout='constrained')
        im=plt.imshow(np.log10(FT2d),**imkwargs,cmap='magma',\
                    extent=[-self.kmax_prime,self.kmax_prime,0,self.wmax_prime],clim=(-18,-13))
        # cbar=plt.colorbar(im)#,orientation='horizontal',shrink=0.6,location='top')
        ax2d.set_xlabel(self.klabel,**tnrfont)
        ax2d.set_ylabel(self.wlabel,**tnrfont)
        # plot solutions ::
        ax2d.plot([0,self.wmax_prime],[0,self.wmax_prime],linestyle='--',color='k')
        ax2d.plot([0,-self.wmax_prime],[0,self.wmax_prime],linestyle='--',color='k')
        # set limits ::
        ax2d.set_xlim(-15,15)
        ax2d.set_ylim(0,self.wmax_prime)
        fig2d.savefig(self.home+'/twoStream_FT2d_Bz.png',bbox_inches='tight')
        plt.show()

if __name__=='__main__':

    home = os.getcwd()

    # sim_loc = getSimulation('/home/space/phrmsf/Documents/EPOCH/epoch-4.17.16/epoch1d/old/old_sims/twoStream')
    sim_loc = getSimulation('/home/space/phrmsf/Documents/EPOCH/5_devel/epoch1d/2stream')

    TSS = tS_Simulation(home,sim_loc)
    TSS.FT1d()
    TSS.FT2d()
