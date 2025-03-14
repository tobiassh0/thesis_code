
from func_load import *


def getKernelDoppler(sims,FT2darr,normspecies,wkmax=[10,25],logthresh=1.8,kernel='custom',dwdkrange=(-0.5,0.5),labels=[None],\
                            theta=86.3,plot=False,home=None,theory_dop=True,kernel_dop=False,dopangle=None):
    """
    Uses the kernel method (see list_new.py and this link [https://pyimagesearch.com/2021/05/12/image-gradients-with-opencv-sobel-and-scharr/]) 
    to find the gradients within the image. Uses a custom/modified Scharr kernel which is more sensitive to 
    -45deg gradients. Plots the histogram of these gradients as well as the converted grad to vA atop the 2d 
    FFT.
        params in
            sims:			list of sims to analyse and extract each gradient from
            FT2darr:		list of 2d FFTs corresponding to the sims in list sims 
            normspecies:the normalisation species used for w,k space
            wkmax:		the limits of the 2d FFT [wmax,kmax] which to plot, in units of normspecies normalisation 
            logthresh:	the log value of the thresh which to apply to the 2d FFT 
            kernel:		the kernel to use (sobel, scharr, custom)
            dwdkrange:	the range over dw/dk (grad) space which to calculate the histogram (modify this for accurate working)
            theta:		redundant parameter, but useful if you want to plot the theoretical doppler as well
            labels:		label with which to annotate the threshed 2d FFTs
        params out
            Plots and saves the histograms of gradients (dw/dk) and 2d FFT with doppler shifted harmonics for the whole 
            array of sims provided. 
            dsvarr:		an array of all doppler velocities in units of vA (per sim)
    """	
    # error check, see if using theory or kernel method
    if theory_dop:
        kernel_dop!=theory_dop # force False for kernel
    if kernel_dop:
        theory_dop!=kernel_dop # force False for theory

    if not home: home = os.getcwd()
    # load sim & FT2d
    l = len(sims)
    if l == 1: l+=1
    if plot:
        fig1,ax1=plt.subplots(nrows=l,figsize=(6,4*len(sims)))
        fig2,ax2=plt.subplots(nrows=3,ncols=l//3,figsize=(9,6),sharex=True)
        fig2.subplots_adjust(hspace=0.1,wspace=0.1)
        ax2=ax2.ravel()
    dsvarr=[] # 1d array of gradient in units of m/s per sim
    for i in range(len(sims)):
        ## calc gradients in image
        # setup
        sim_loc = getSimulation(sims[i])	
        print(sim_loc)
        d0 = sdfread(0)
        times = read_pkl('times')
        vA = getAlfvenVel(d0)
        print(vA/const.c)
        klim = 0.5*2*const.PI/getdxyz(d0)
        wlim = 0.5*2*const.PI/getdt(times)
        wnorm = getCyclotronFreq(d0,normspecies)
        knorm = wnorm/vA
        # freq resolutions (no factor half)
        dk = 2*const.PI/getGridlen(d0)
        dw = 2*const.PI/times[-1]
    
        # cut FT2d into size needed
        FT2d = FT2darr[i]
        (nw,nk) = FT2d.shape
        wmax = wkmax[0]*wnorm ; kmax = wkmax[1]*knorm
        FT2d = np.log10(FT2d[:int(nw*wmax/wlim),:int(nk*kmax/klim)])
        (nw,nk) = FT2d.shape

        # threshold FT2d
        tFT2d= FT2d.copy() # copy arr
        tFT2d[FT2d < logthresh] = 0
        ttFT2d = tFT2d.copy()
        ttFT2d[FT2d > logthresh] = 1
        FT2d = ttFT2d # overwrite old FT2d 
        del ttFT2d # delete temp arr
        # plt.imshow(FT2d,aspect='auto',origin='lower',extent=[0,wkmax[0],0,wkmax[1]])
        
        # check if using theory or kernel estimation for Doppler angle
        if theory_dop:
            """ theory doppler line """
            Ep = const.qe*14.68e6
            up = np.sqrt(2*Ep/(const.me_to_mp*const.me))
            uperp = 0.9*vA
            pitch_angle = np.arcsin(uperp/up)
            dShiftvel = (up/vA)*np.cos(89*const.PI/180)*np.cos(pitch_angle)
        elif kernel_dop:
            """ Kernel gradient estimation (empirical) """
            # Kernel gradient map
            _,kGangle = Kernel(FT2d,kernel=kernel) # list_new func : kernel = 'scharr' or 'sobel'

            # gradients as angles
            kGangle = kGangle[1:-1,1:-1]# remove abberations around edge
            # plt.imshow(kGangle,extent=[0,wkmax[1],0,wkmax[0]],**kwargs)	; plt.show()	
            kGangle = kGangle.flatten()

            # convert to all negative angles (easier to calc real gradient)
            for g in range(len(kGangle)):
                if kGangle[g] > 0:
                    kGangle[g]-=const.PI #rad

            # convert grad from angle to velocity (convert +-inf to nan)
            dwdk = np.nan_to_num(np.tan(kGangle),posinf=np.nan,neginf=np.nan)
            # remove nan & zero values
            dwdk = dwdk[~np.isnan(dwdk)]
            dwdk = dwdk[dwdk!=0]

            # normalise & thresh
            dw_dk = dwdk * (dw/dk)/vA # normalised to vA
            #thresh = (np.abs(dw_dk) < dwdkrange[1])
            #dw_dk = dw_dk[thresh]
            print(kernel+' kernel mean :: ',np.mean(dw_dk))
            print(kernel+' kernel medi :: ',np.median(dw_dk))

            # calc histogram
            #ax1[i].hist(counts,bins=bins) dw_dk
            counts,bins=np.histogram(dw_dk,bins=1000,density=True,range=dwdkrange)
            dShiftvel = bins[np.argmax(counts)] # doppler shift velocity, in units of vA
            dsvarr.append([dsv*vA,dsv,vA])
            print(kernel+' kernel max :: ', dsv)
        else:
            """ other method - given Doppler angle """
            if dopangle==None:
                print('# ERROR # : No Doppler angle given')
                raise SystemError
            else:
                dShiftvel=((dw/dk)*(1/np.tan(dopangle)))/vA #angle to velocity equivalent (normalised to vA)

        print('Doppler velocity [vA] :: ',dShiftvel)

        # plotting script
        if plot:
            # # plot hist
            # counts,bins,_=ax1[i].hist(dw_dk,bins=1000,density=True,range=dwdkrange) # np.log10
            # ax1[i].set_ylabel('Normalised count',**tnrfont)
            # #dsva.append([dsv,vA])

            # doppler shifted line 
            kx = np.linspace(0,20,100)*knorm
            for j in range(0,int(wmax/wnorm),1):
                w = (j*wnorm)*np.ones(len(kx))
                # plot Doppler shifted lines for each ICE harmonic
                ww = w - (abs(dShiftvel)*vA)*kx
                ax2[i].plot(kx/knorm,ww/wnorm,color='darkgrey',linestyle='--',zorder=1)

            # plot FT2d
            tFT2d = np.ma.masked_where(tFT2d == 0, tFT2d)
            im = ax2[i].imshow((tFT2d),**imkwargs,extent=[0,wkmax[1],0,wkmax[0]],cmap='magma',zorder=2)
            del tFT2d

            # # vA line
            # ax2[i].plot([0,10],[0,10],color='k',linestyle='--')

            # annotate concentrations
            ax2[i].annotate(hlabels[i],xy=(0.1,0.98),xycoords='axes fraction',color='k',ha='left',va='top')

            # ax limits
            ax2[i].set_xlim(0,20) # reduce plotting limits
            ax2[i].set_ylim(0,10) # " " 

            # set-ticks and ax formatting
            if i >= 2*len(sims)//3:
                ax2[i].set_xticks([0,5,10,15]) ; ax2[i].set_xticklabels([0,5,10,15])
            else:
                ax2[i].set_xticks([]) ; ax2[i].set_xticklabels([])
            if i%3==0:
                ax2[i].set_yticks([0,2,4,6,8,10]) ; ax2[i].set_yticklabels([0,2,4,6,8,10])
            else:
                ax2[i].set_yticks([]) ; ax2[i].set_yticklabels([])

    # total plot formatting
    if plot:
        ## format edge axes
        fig2.supylabel(r'$\omega/$'+getOmegaLabel(normspecies),**tnrfont,x=0.05)
        fig2.supxlabel(r'$kv_A/\Omega_p$',**tnrfont,y=-0.03)

        # ax2[len(ax2)//2].set_ylabel(r'$\omega/$'+getOmegaLabel(normspecies),**tnrfont)
        # ax2[len(ax2)//2].set_yticks([0,2,4,6,8,10]) ; ax2[len(ax2)//2].set_yticklabels([0,2,4,6,8,10])
        ax2[-1].set_xticks([0,5,10,15,20]) ; ax2[-1].set_xticklabels([0,5,10,15,20])
        ax1[-1].set_xlabel(r'$d\omega/dk$'+'  '+r'$[v_A]$',**tnrfont)
        # fig1.savefig('dw_dk_'+kernel+'_grad.png',bbox_inches='tight')
        # ax2[-1].set_xlabel(r'$kv_A/\Omega_p$',**tnrfont)

        # add colorbar to ax2 (FT2d)
        p00 = ax2[0].get_position().get_points().flatten()
        p22 = ax2[-1].get_position().get_points().flatten()
        ax2_cbar = fig2.add_axes([p00[0], 0.97, p22[2]-p00[0], 0.02]) # [left bottom width height]
        cbar = plt.colorbar(im, cax=ax2_cbar, orientation='horizontal')
        fig2.savefig(home+'/FT_2d_doppler_all.png',bbox_inches='tight')
        # plt.show()
    return None # np.array(dsvarr)

def PLOTDOPPLER(hlabels,dsvarr,pchange=True):
    """
    Plotting the (1) Alf speed vA (2) Kernel extracted Doppler velocity and (3) the Doppler velocity
    normalised to the 0% He3 concentration Doppler shift (i.e. vA normalisation removed)
        IN:
            hlabels : the He3 concentration as a % (e.g. 0, 5, 10, 15 ...)
            dsvarr 	: an array of the Doppler extracted velocities (dsv) and vA for each sim
            pchange : boolean whether to plot percentage change (all on one plot) or with nrows
        OUT:
            Saves the figure 
    """
    hlabels = hlabels/100 # change to actual value not %
    dsv_vA,dsv,vA = dsvarr[:,0], dsvarr[:,1], dsvarr[:,2]
    if pchange:
        fig,ax=plt.subplots(figsize=(8,8/const.g_ratio))
        # percentage change
        ax.scatter(hlabels,vA/vA[0]-1,color='k',marker='o')
        ax.scatter(hlabels,abs(dsv/dsv[0])-1,color='k',marker='x')
        ax.scatter(hlabels,abs(dsv_vA)/abs(dsv_vA[0])-1,color='k',marker='s')
        # fit linear curves
            # va
        params, params_err = ODR_fit(x=hlabels,y=vA/vA[0]-1)
        ax.plot(hlabels,func_linear(params,hlabels),color='r',linestyle='--')
        top = (params[0] + params_err[0], params[1] + params_err[1])
        bottom = (params[0] - params_err[0], params[1] - params_err[1])
        ax.fill_between(hlabels,hlabels*top[0]+top[1],hlabels*bottom[0]+bottom[1],color='r',alpha=1/3)
            # vdop
        params, params_err = ODR_fit(x=hlabels,y=abs(dsv/dsv[0])-1)
        ax.plot(hlabels,func_linear(params,hlabels),color='r',linestyle='--')
        top = (params[0] + params_err[0], params[1] + params_err[1])
        bottom = (params[0] - params_err[0], params[1] - params_err[1])
        ax.fill_between(hlabels,hlabels*top[0]+top[1],hlabels*bottom[0]+bottom[1],color='r',alpha=1/3)
            # vdop %
        params, params_err = ODR_fit(x=hlabels,y=abs(dsv_vA)/abs(dsv_vA[0])-1)
        ax.plot(hlabels,func_linear(params,hlabels),color='r',linestyle='--')
        top = (params[0] + params_err[0], params[1] + params_err[1])
        bottom = (params[0] - params_err[0], params[1] - params_err[1])
        ax.fill_between(hlabels,hlabels*top[0]+top[1],hlabels*bottom[0]+bottom[1],color='r',alpha=1/3)
        # formatting
        ax.set_xlim(-0.05,0.5)
        ax.set_ylim(-0.15,0.15)
        ax.set_xlabel(r'$\xi_{He3}$',**tnrfont)
        ax.set_ylabel('Fractional change',**tnrfont)
        ax.legend(labels=[r'$v_A$',r'$v_{dop}/v_A$',r'$v_{dop}$'],loc='best',fontsize=18)
        fig.savefig('vA_vdop_fits_pchange.png',bbox_inches='tight')
    else:
        fig,axs=plt.subplots(figsize=(4,8),nrows=3,sharex=True)
        fig.subplots_adjust(hspace=0.1)
        ##	vA
        axs[0].scatter(hlabels,vA/const.c,color='k')
        axs[0].set_ylabel(r'$v_A/c$',**tnrfont)
        axs[0].set_ylim(0.0265,0.031)
        # fit linear curve
        params, params_err = ODR_fit(x=hlabels,y=vA/const.c)
        axs[0].plot(hlabels,func_linear(params,hlabels),color='r',linestyle='--')
        top = (params[0] + params_err[0], params[1] + params_err[1])
        bottom = (params[0] - params_err[0], params[1] - params_err[1])
        axs[0].fill_between(hlabels,hlabels*top[0]+top[1],hlabels*bottom[0]+bottom[1],color='r',alpha=1/3)
        # pearsons cross cor
        r = np.corrcoef(hlabels,vA/const.c)[0,1]
        axs[0].annotate(r'$r={}$'.format(np.around(r,3)),xy=(0.03,0.85),xycoords='axes fraction',ha='left',va='bottom',fontsize=18)
        ## v_dop
        axs[1].scatter(hlabels,abs(dsv),color='k')
        axs[1].set_ylabel(r'$|v_{dop}/v_A|$',**tnrfont)
        axs[1].set_ylim(0.118,0.138)
        # fit linear curve
        params, params_err = ODR_fit(x=hlabels,y=abs(dsv))
        axs[1].plot(hlabels,func_linear(params,hlabels),color='r',linestyle='--')
        top = (params[0] + params_err[0], params[1] + params_err[1])
        bottom = (params[0] - params_err[0], params[1] - params_err[1])
        axs[1].fill_between(hlabels,hlabels*top[0]+top[1],hlabels*bottom[0]+bottom[1],color='r',alpha=1/3)
        # pearsons cross cor
        r = np.corrcoef(hlabels,abs(dsv))[0,1]
        axs[1].annotate(r'$r={}$'.format(np.around(r,3)),xy=(0.97,0.85),xycoords='axes fraction',ha='right',va='bottom',fontsize=18)
        ## v_dop %
        axs[2].scatter(hlabels,abs(dsv_vA)/abs(dsv_vA[0]),color='k')
        axs[2].set_ylabel(r'$|v_{dop}/v_{dop}^{(0)}|$',**tnrfont)
        axs[2].set_xlabel(r'$\xi_{He3}$',**tnrfont)
        axs[2].set_ylim(0.95,1.05)
        # fit linear curve
        params, params_err = ODR_fit(x=hlabels,y=abs(dsv_vA)/abs(dsv_vA[0]))
        axs[2].plot(hlabels,func_linear(params,hlabels),color='r',linestyle='--')
        top = (params[0] + params_err[0], params[1] + params_err[1])
        bottom = (params[0] - params_err[0], params[1] - params_err[1])
        axs[2].fill_between(hlabels,hlabels*top[0]+top[1],hlabels*bottom[0]+bottom[1],color='r',alpha=1/3)
        # pearsons cross cor
        r = np.corrcoef(hlabels,abs(dsv_vA)/abs(dsv_vA[0]))[0,1]
        axs[2].annotate(r'$r={}$'.format(np.around(r,3)),xy=(0.97,0.85),xycoords='axes fraction',ha='right',va='bottom',fontsize=18)
        axs[2].set_xlim(-0.05,0.5)
        fig.savefig('vA_vdop_fits.png',bbox_inches='tight')
    ## 
    plt.show()
    return None

def doppler_onesim(sim='/storage/space2/phrmsf/lowres_D_He3/0_38_p_90'):
    """
    Example function of plotting the Doppler lines against a FT2d for fixed minority species 
    parameters. 
    """
    # load sim & FT2d 
    sim_loc = getSimulation(sim)
    FT2d = read_pkl('FT_2d_Magnetic_Field_Bz')
    d0 = sdfread(0)
    times = read_pkl('times')
    vA = getAlfvenVel(d0)
    print(vA)
    Ep = 14.68*1e6*const.qe
    vp = (2*Ep/getMass('Protons'))**0.5 # proton velocity
    klim = 0.5*2*const.PI/getdxyz(d0) 
    wlim = 0.5*2*const.PI/getdt(times)
    wcp = getCyclotronFreq(d0,'Protons')
    wcD = getCyclotronFreq(d0,'Deuterons')
    wnorm = wcp
    knorm = wnorm/vA
    (nw,nk) = (FT2d.shape)
    dk = 2*const.PI/getGridlen(d0) # no factor half (see e-notes 24/10/23)
    dw = 2*const.PI/times[-1]
    print(dw,dk)
    # figure setup
    fig,ax=plt.subplots(ncols=3,figsize=(6*3,4))
    
    ## calc gradients in image
    # cut FT2d into size needed
    (nw,nk) = FT2d.shape
    kmax = 20*knorm ; wmax = 10*wnorm
    kmin = 0*knorm; wmin = 0#*wnorm
    FT2d = np.log10(FT2d[:int(nw*wmax/wlim),int(nk*kmin/klim):int(nk*kmax/klim)])
    (nw,nk) = FT2d.shape
    
    # threshold FT2d
    thresh = 1.8
    tFT2d = FT2d.copy()
    #tFT2d[FT2d > 1.6] = 0
    tFT2d[FT2d < thresh] = 0
    
    # Kernel gradient map
    kernel = 'custom'
    _,kGangle = Kernel(tFT2d,kernel=kernel) # scharr, sobel or custom
    # gradients as angles
    kGangle = kGangle[1:-1,1:-1]# remove abberations around edge
    im = ax[0].imshow(kGangle*180/const.PI,**kwargs,extent=[0,kmax/knorm,0,wmax/wnorm],cmap='Accent')
    plt.colorbar(im)
    print((dw/dk)/vA)
    # remove zero values (dont want to plot them in hist)
    dwdk = np.nan_to_num(np.tan(kGangle),posinf=np.nan,neginf=np.nan) # remove inf and -inf values
    dw_dk = dwdk * (dw/dk)/vA # normalise
    kGangle = kGangle.flatten()
    # convert to all negative angles (easier to calc real gradient)
    for i in range(len(kGangle)):
        if kGangle[i] > 0:
            kGangle[i]-=const.PI #rad
    # remove inf values
    dwdk = np.nan_to_num(np.tan(kGangle),posinf=np.nan,neginf=np.nan) # remove inf and -inf values
    dwdk = dwdk[~np.isnan(dwdk)]
    # remove zero values (dont want to plot them in hist)
    dwdk = dwdk[dwdk!=0]
    dw_dk = dwdk * (dw/dk)/vA # normalise
    thresh = (np.abs(dw_dk) < 2.0) & (np.abs(dw_dk) > 0.001)
    dw_dk = dw_dk[thresh]
    print(kernel+' kernel mean :: ',np.mean(dw_dk))
    print(kernel+' kernel medi :: ',np.median(dw_dk))
    # plot hist
    counts,bins,_=ax[1].hist(dw_dk,bins=1000,density=True,range=(-1,1)) # np.log10
    dsv = bins[np.argmax(counts)] # doppler shift velocity in units of vA
    print(kernel+' kernel max :: ', dsv)
    ax[1].set_ylabel('Normalised count',**tnrfont)
    ax[1].set_xlabel(r'$d\omega/dk$'+'  '+r'$[v_A]$',**tnrfont)
    # plot kde
    kde = stats.gaussian_kde(dw_dk)
    xx = np.linspace(-1,1,1000)
    ax[1].plot(xx,kde(xx),color='r')
    #fig.savefig('dw_dk_'+kernel+'_grad.png',bbox_inches='tight')
    #plt.show()
    
    # plot FT2d
    #fig,ax=plt.subplots(figsize=(8,6))
    ax[2].imshow((tFT2d[1:,1:]),**kwargs,extent=[0,kmax/knorm,0,wmax/wnorm])
    kx = np.linspace(0,20,100)*knorm
    theta = 86.3 # deg
    kperp = np.sin(theta*const.PI/180)
    kpara = np.cos(theta*const.PI/180)
    uperp = 0.9; upara = 6.076 
    dsth = -(kperp*uperp + kpara*upara)
    print(dsv,(dsth+1))
    # doppler shifted line 
    for i in range(0,int(wmax/wnorm),1):
        w = wnorm*np.ones(len(kx))*i
        # empirical
        ww = w + (dsv*vA)*kx
        ax[2].plot(kx/knorm,ww/wnorm,color='white',linestyle='--')
    #	# theory
    #	tww = w + ((dsth+1)*vA)*kx
    #	ax.plot(kx/knorm,tww/wnorm,color='white',linestyle='-.')
    ax[2].set_xlim(0,20)
    ax[2].set_ylim(0,20)
    ax[2].set_ylabel(r'$\omega/\Omega_p$',**tnrfont)
    ax[2].set_xlabel(r'$kv_A/\Omega_p$',**tnrfont)
    ax[2].plot([0,10],[0,10],color='white',linestyle='--') # vA line
    
    #fig.savefig('FT_2d_doppler_th.png',bbox_inches='tight')
    plt.show()

if __name__ == '__main__':
    import scipy
    import matplotlib

    # # D-He3
    quantity='Magnetic_Field_Bz'
    home = '/storage/space2/phrmsf/lowres_D_He3/'
    sims = np.sort([i for i in os.listdir(home) if 'p_90' in i])
    hlabels = np.array([int(i[2:4]) for i in sims])
    # collect all FT2d (needed for Doppler shifts)
    FT2darr = []
    for sim in sims:
        _=getSimulation(home+sim)
        FT2darr.append(read_pkl('FT_2d_'+quantity))
    os.chdir(home)

    # kernel
    sims = [home+sim for sim in sims]
    dsvarr = getKernelDoppler(sims,FT2darr,labels=hlabels,normspecies='Protons',plot=True)
    # PLOTDOPPLER(hlabels,dsvarr)	
        