
from func_load import *
if __name__=='__main__':
    home='/storage/space2/phrmsf/lowres_D_He3/'
    sims=np.sort([i for i in os.listdir(home) if 'p_90' in i])
    xiHe3=[int(i.split('_')[1]) for i in sims]
    omegas_match=np.zeros(len(sims))
    fig,axs=plt.subplots(figsize=(15,8),nrows=3,ncols=3,sharex=True,sharey=True,layout='constrained')
    fig_sep,ax_sep=plt.subplots(figsize=(6,6),layout='constrained')
    ax=axs.ravel()
    i=0
    for sim in sims:
        simloc=getSimulation(home+sim)
        # load 0th file (cold plasma disp)
        d0=sdfread(0)
        times=read_pkl('times')
        # dispersion limits & norm
        knyq=0.5*2*const.PI/getdxyz(d0)
        wnyq=0.5*2*const.PI/getdt(times)
        dw=2*wnyq/len(times)
        wcp=getCyclotronFreq(d0,'Protons')
        vA=getAlfvenVel(d0)
        # load FT2d data
        FT2d=read_pkl('FT_2d_Magnetic_Field_Bz')
        # ax[i].imshow(np.log10(FT2d),**imkwargs,cmap='magma',clim=(-4,6),extent=[0,knyq*vA/wcp,0,wnyq/wcp])
        # cold plasma dispersion
        omegas=wcp*np.linspace(0,20,int(1e6))
        k1,k2,k3=coldplasmadispersion(d0,omegas,theta=None)
        # # solution overlap plot
        # ax[i].plot(np.imag(k1)*vA/wcp,omegas/wcp,'w.')
        # ax[i].plot(np.real(k2)*vA/wcp,omegas/wcp,'k,')        
        # find where solutions overlap (greater than 10wcp to make easier)
        thresh=omegas/wcp>10
        index=np.argmin(np.abs(np.imag(k1[thresh])-np.real(k2[thresh])))
        # load & plot power spectra
        log10_power=read_pkl('log10_power')
        omegas_power=read_pkl('omegas_power')
        ax[i].plot(omegas_power/wcp,log10_power-np.log10(dw),color='k')
        # plot vertical line of evanescent wave intersection
        omega_match=omegas[thresh][index]
        omegas_match[i]=omega_match/wcp
        ax[i].annotate(str(xiHe3[i]),xy=(0.05,0.95),xycoords='axes fraction',va='top',ha='left')
        for j in range(0,21):
            ax[i].axvline(j,color='darkgrey',linestyle='--')
        ax[i].axvline(omega_match/wcp,linestyle='--',color='r')
        i+=1
    ax_sep.plot(np.array(xiHe3)/100,omegas_match,'ro-')
    ax_sep.plot([0,0.45],[omegas_match[0],omegas_match[-1]],'k-')
    ax[-1].set_xlim(0,20)
    ax_sep.set_ylabel(r'$\omega_{ev.}/\Omega_p$',**tnrfont)
    ax_sep.set_xlabel(r'$\xi_{He3}$',**tnrfont)
    fig.supylabel('PSD',**tnrfont) ; fig.supxlabel(r'$\omega/\Omega_p$',**tnrfont)
    fig.savefig(home+'evanescent_match_PSD.png')
    fig_sep.savefig(home+'evanescent_match_freq.png')
    # plt.show()