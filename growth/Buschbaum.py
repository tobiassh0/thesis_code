
from func_load import *
import numpy as np
import matplotlib.pyplot as plt
import my_constants as const

class Buschbaum(object):
    def __init__(self,home,species,_B0=2.1,_ne=1e19,_xi3=2e-3):
        # constants
        self.home = home
        e0 = const.e0
        qe = const.qe
        self.B0 = _B0
        self.ne = _ne
        self.xi3 = _xi3
        # charge and masses
        self.Z1 = getChargeNum(species[0])
        self.Z2 = getChargeNum(species[1])
        self.Z3 = getChargeNum(species[2])
        self.Me = getMass('Electrons')
        self.M1 = getMass(species[0])
        self.M2 = getMass(species[1])
        self.M3 = getMass(species[2])
        # frequencies
        self.wp = np.sqrt((self.ne*(qe)**2)/(self.Me*e0))
        self.wb = qe*self.B0/self.Me
        self.wb1 = self.Z1*qe*self.B0/self.M1
        self.wb2 = self.Z2*qe*self.B0/self.M2

    def setupBB(self,_xi2):
        # (relative) concentrations
        if np.all(_xi2==None):
            xi2 = np.array([0.1,0.2,0.5,0.75,0.85,0.95])
        else:
            xi2 = np.array(_xi2) # np.linspace(0,1,100)
        xi1 = (1/self.Z1)*(1-self.Z2*xi2-self.Z3*self.xi3)
        rel_xi1 = xi1/(xi1+xi2) # relative concentration of ions
        rel_xi2 = xi2/(xi1+xi2)
        if np.all(rel_xi1+rel_xi2 != 1):
            print(rel_xi1+rel_xi2) # array of ones, len (xi2)
            print('# ERROR # :: Index [{}] of relative concentrations don\'t add to 1'.format(np.where(rel_xi1+rel_xi2 != 1)[0][0]))
            raise SystemError
        # reformat and update values
        self.x1 = rel_xi1
        self.x2 = rel_xi2
        del rel_xi1, rel_xi2 
        self.f1 = self.Me/self.M1
        self.f2 = self.Me/self.M2
        pass

    def heatmap(self,_xi2=None,omega_max=4.0):
        # get concentrations and mass ratios
        self.setupBB(_xi2)
        # frequencies
        omega = self.wb2*np.linspace(0,omega_max,1000)
        # setup figure
        fig,ax=plt.subplots(figsize=(8,8))
        # meshgrid for heatmap (w, xi2, krel)
        X2, OMEGA = np.meshgrid(self.x2,omega)
        X1 = 1-X2
        twoion = OMEGA**2 - ((X1*self.f2+X2*self.f1)/(X1*self.f1+X2*self.f2))*self.wb1*self.wb2
        # numerator
        A = self.wp**2 * (OMEGA**2 - (X2*self.wb1 + X1*self.wb2)**2)
        B = (OMEGA**2 - self.wb*(X1*self.wb1 + X2*self.wb2))*twoion
        num = A - B
        # denominator
        C = (OMEGA**2 - self.wb**2)*(OMEGA**2 - self.wb1**2)*(OMEGA**2 - self.wb2**2)
        D = (self.wp**2 * (OMEGA**2 - self.wb*(X1*self.wb1 + X2*self.wb2)))*twoion
        den = C - D
        # wavenumber, normalised and relative
        krel = np.sqrt(1+(self.wp**2)*(num/den))
        # plot
        im=ax.imshow(np.log10(krel),**imkwargs,cmap='magma',\
                    extent=[self.x2[0],self.x2[-1],omega[0]/self.wb2,omega[-1]/self.wb2])
        # thresh = (krel>0) & (krel<3500)
        # im=ax.imshow(np.log10(krel*thresh),**imkwargs,cmap='magma',\
        #     extent=[self.x2[0],self.x2[-1],omega[0]/self.wb2,omega[-1]/self.wb2])            
        cbar=plt.colorbar(im)
        cbar.set_label(label=r'$\log_{10}(k_r)$',fontsize=18)
        ax.set_ylabel(r'$\omega/\Omega_2$',**tnrfont)
        ax.set_xlabel(r'$\xi_2$',**tnrfont)
        fig.savefig(home+'/buschbaum_heatmap.png')
        pass

    def line_plot(self,_xi2=None,omega_max=4.0):
        # setup concentrations and mass ratios
        self.setupBB(_xi2)
        # frequencies
        omega = self.wb2*np.linspace(0,omega_max,1000)
        # setup figure
        fig,ax=plt.subplots(figsize=(8,8))
        # line plot
        for i in range(len(self.x1)):
            _x1 = self.x1[i] ; _x2 = self.x2[i]
            twoion = omega**2 - ((_x1*self.f2+_x2*self.f1)/(_x1*self.f1+_x2*self.f2))*self.wb1*self.wb2
            # numerator
            A = self.wp**2 * (omega**2 - (_x2*self.wb1 + _x1*self.wb2)**2)
            B = (omega**2 - self.wb*(_x1*self.wb1 + _x2*self.wb2))*twoion
            num = A - B
            # denominator
            C = (omega**2 - self.wb**2)*(omega**2 - self.wb1**2)*(omega**2 - self.wb2**2)
            D = (self.wp**2 * (omega**2 - self.wb*(_x1*self.wb1 + _x2*self.wb2)))*twoion
            den = C - D
            # wavenumber, normalised and relative
            krel = np.sqrt(1+(self.wp**2)*(num/den))
            # frequency of divergent wavenumber
            omega_div = omega[np.argmax(krel)]
            # frequency of wavenumber at turnover (0 grad)
            omega_turn=omega[np.argmax(np.gradient(krel,(omega[-1]-omega[0])/len(omega)))]
            ax.axvline(omega_div/self.wb2,linestyle='--',color='b')
            ax.axvline(omega_turn/self.wb2,linestyle='--',color='r')
            # plot
            thresh = [True]*len(krel) # krel>0
            ax.plot(omega/self.wb2,krel,color='k',label=r"{:.2f}%".format(100*_x2))
            # ax.plot(omega[thresh]/wb2,krel[thresh],',',color='k',label=r"{:.2f}%".format(100*x2[i]))
        # ax.legend(loc='best')
        ax.set_xlabel(r'$\omega/\Omega_2$',**tnrfont)
        ax.set_ylabel(r'$k_r$',**tnrfont)
        # fig.savefig(self.home+'/buschbaum_line.png')
        plt.show()
        pass

# ================================ #

if __name__=='__main__':
    # home dir
    home = os.getcwd()
    bb = Buschbaum(home,['Deuterons','Tritons','Alphas'])
    bb.line_plot(_xi2=[0.5])
    # bb.heatmap()
    