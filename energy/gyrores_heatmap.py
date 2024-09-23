
from func_load import *

if __name__=='__main__':
    home = '/storage/space2/phrmsf/lowres_D_He3/'
    sims = [i for i in os.listdir(home) if 'p_90' in i]
    mD = getMass('Deuterons') ; mHe3 = getMass('He3')
    qqD = (const.qe)**2 ; qqHe3 = (2*const.qe)**2
    for sim in sims:
        print(sim)
        simloc=getSimulation(home+sim)
        times=read_pkl('times')
        T=times[-1]
        d0=sdfread(0)
        xiD,xiHe3,_=getConcentrationRatios(d0)
        print(xiD,xiHe3)
        tcp=2*const.PI/getCyclotronFreq(d0,'Protons')
        L=getGridlen(d0)
        all_species=getIonSpecies(d0)
        fig,axs=plt.subplots(figsize=(6,12),nrows=3)
        axs=axs.ravel()
        if 'He3' in all_species:
            # load Energies in (X,T)
            uD_xt=read_pkl('Deuterons_KEdensmatrix')
            uHe3_xt=read_pkl('He3_KEdensmatrix')
            # change in energy components
            duD_xt = uD_xt - np.mean(uD_xt[0:10,:])
            duHe3_xt = uHe3_xt - np.mean(uHe3_xt[0:10,:])
            axs[0].imshow(duD_xt,cmap='magma',**imkwargs,extent=[0,L,0,T/tcp])
            axs[1].imshow(duHe3_xt,cmap='magma',**imkwargs,extent=[0,L,0,T/tcp])
            im = axs[2].imshow((duD_xt/duHe3_xt)/((mHe3*qqD*xiD)/(mD*qqHe3*xiHe3)),cmap='bwr',**imkwargs,extent=[0,L,0,T/tcp],clim=(0,2))
            cbar=plt.colorbar(im)
            fig.supylabel('Time '+r'$[\tau_{cp}]$',**tnrfont)
            fig.supxlabel(r'$x$'+' [m]',**tnrfont)
            fig.savefig(home+sim+'/gyrores_ratio.png',bbox_inches='tight')
            plt.show()
        else:
            pass