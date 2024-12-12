
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

    def heatmap(self,_xi2=None,omega_max=None):
        # if haven't set frequency limit
        if not omega_max:
            omega_max = float(self.M2/self.M1)
        # get concentrations and mass ratios
        self.setupBB(_xi2)
        # frequencies
        omega = self.wb2*np.linspace(0,omega_max,1000)
        # setup figure
        fig,ax=plt.subplots(figsize=(6,6),layout='constrained')
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

    def steady_stateBB(self,omega,krel,omega_div,diffkgrad=0.01):
        """
            Finds the value of the frequency (less than the divergent frequency)
            which is considered the turnover/divergence from the steady-state.
            In:
                self        : class object 
                omega       : np.array of frequencies (unnormalised) used to calculate krel
                krel        : np.array of wavenumber solutions, normalised by k0=w/c
                omega_div   : value of frequency where krel solutions diverge
                diffkgrad   : float value percentage (%/100) of gradient in krel which to consider 
                                loss of steady state
            Out:
                omega_turn  : float value of turnover frequency (unnormalised) 
        """
        # get gradient of wavenumber
        kgrad = np.gradient(krel)
        # frequency of turnover < diverging
        thresh = omega<omega_div 
        # 1% floor value (turnover from steady-state)
        index = np.nanargmin(np.abs(kgrad[thresh]-diffkgrad))
        return (omega[thresh])[index]

    def line_plot(self,_xi2=None,omega_max=None):
        # if haven't set frequency limit
        if not omega_max:
            omega_max = float(self.M2/self.M1)
        # setup concentrations and mass ratios
        self.setupBB(_xi2)
        # frequencies
        omega = self.wb2*np.linspace(0,omega_max,1000)
        # setup figure
        fig,ax=plt.subplots(figsize=(6,6),layout='constrained')
        fig_ra,ax_ra=plt.subplots(figsize=(6,6),layout='constrained')
        res_amp=np.zeros(len(self.x1))
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
            omega_turn = self.steady_stateBB(omega,krel,omega_div) # omega[np.argmax(np.gradient(krel,(omega[-1]-omega[0])/len(omega)))]
            # ax.axvline(omega_div/self.wb2,linestyle='--',color='b') # divergent frequency
            # ax.axvline(omega_turn/self.wb2,linestyle='--',color='r') # turnover frequency (steady-state)
            res_amp[i]=(omega_div-omega_turn)/self.wb2
        #     # plot
        #     thresh = [True]*len(krel) # krel>0
        #     ax.plot(omega/self.wb2,krel,color='k',label=r"{:.2f}%".format(100*_x2))
        #     # ax.plot(omega[thresh]/wb2,krel[thresh],',',color='k',label=r"{:.2f}%".format(100*x2[i]))
        # # line plot
        # ax.set_xlabel(r'$\omega/\Omega_2$',**tnrfont)
        # ax.set_ylabel(r'$k_r$',**tnrfont)
        # fig.savefig(self.home+'/buschbaum_line.png')
        # resonance amplitude
        ax_ra.plot(self.x2,res_amp,'-o',color='k')
        ax_ra.axvline(self.x2[np.argmax(res_amp)],linestyle='--',color='r')
        ax_ra.text(self.x2[np.argmax(res_amp)],0.95*np.max(res_amp), "{:.2f}".format(self.x2[np.argmax(res_amp)]),color='r',\
                    ha='right', va='top', rotation=90,**tnrfont)
        ax_ra.set_ylabel(r'Resonance Amplitude, '+r'$(\omega_\infty-\omega_{\delta})$',**tnrfont)
        ax_ra.set_xlabel(r'Relative Concentration, '+r'$\xi_2/(\xi_1+\xi_2)$',**tnrfont)
        fig_ra.savefig(self.home+'/buschbaum_resamp.png')
        pass

    def maxresamp(self):

        pass

# ================================ #

if __name__=='__main__':
    import csv
    # home dir
    home = os.getcwd()
    f = open('bb-resamp.txt','r')
    maxresamp=[] ; massra=[]
    for x in f:
        a,b=x.strip("\n").split(" ")
        massra.append(float(a)) ; maxresamp.append(float(b))
    fig,ax=plt.subplots(figsize=(8,5),layout='constrained')
    ax.plot(massra[1:],maxresamp[1:],'-o',color='k')
    ax.set_ylabel(r'Concentration of Max. Res. Amp.',**tnrfont)
    ax.set_xlabel(r'Ion Mass Ratio, '+r'$(M_2/M_1)$',**tnrfont)
    ax.set_ylim(0)
    fig.savefig(home+'/buschbaum_maxresamp.png')
    plt.show()
    sys.exit()

    # class var
    bb = Buschbaum(home,['Deuterons','Tritons','Alphas'])
    # plots (line, heatmap & max resonant amplitude)
    bb.line_plot(_xi2=np.linspace(0,0.95,100))
    # bb.heatmap(_xi2=None)
    