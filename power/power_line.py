
from func_load import *
import scipy
import dopplerShift.line_integrate as li


def getPowerLine(FT2d,ynarr,freqs,kmax_prime=100,wmax_prime=20,angle=-3.14/2):
    Nw, Nk = FT2d.shape
    # find power along line
    power_line = np.zeros(len(freqs))
    lsizearr = []
    for i in range(len(ynarr)):
        yn = ynarr[i]
        xlim,ylim = li.LineBoxIntersection(Ys=Nw,Ye=0,Xn=0,Yn=yn,XL=Nk,YL=Nw,\
                                            theta=angle)
        # number of cells within line
        lsizearr.append(np.around(np.sqrt(((xlim[1]-xlim[0]))**2 + ((ylim[1]-ylim[0]))**2))) # No. pixels
        # pixel coordinate points along line
        x = np.linspace(xlim[0],xlim[1],int(lsizearr[-1]))
        y = np.linspace(ylim[0],ylim[1],int(lsizearr[-1]))
        # values along line
        zi = scipy.ndimage.map_coordinates(FT2d,np.vstack((y,x)),order=0)
        # power along line
        power_line[i] = np.sum(zi**2)
    return power_line

def plotPowerLine(home,sim,xi2,wkmax=[10,20],minspec='Protons',maxang=-90,defangle=-85):
    simloc = getSimulation(sim)
    # load FFT2d and plot
    FT2d = read_pkl('FT_2d_Magnetic_Field_Bz')
    Nw, Nk = FT2d.shape
    times = read_pkl('times')
    klim, wlim = getDispersionlimits(sim)
    kmin, kmax = klim
    wmin, wmax = wlim
    # frequency and Alfven speed
    wcmin = getCyclotronFreq(sdfread(0),minspec)
    vA = getAlfvenVel(sdfread(0))
    # normalised cell widths
    dk = (2*const.PI/getGridlen(sdfread(0)))
    dw = (2*const.PI/times[-1])
    # extents = [0,kmax*vA/wcmin,0,wmax/wcmin]
    # thresh FT2d
    wmax_prime, kmax_prime = wkmax
    extents = [0,kmax_prime,0,wmax_prime]
    # threshold FT2d
    kthresh = int(FT2d.shape[1]*extents[1]*wcmin/vA/kmax)
    wthresh = int(FT2d.shape[0]*extents[-1]*wcmin/wmax)
    tFT2d = FT2d[:wthresh,:kthresh]
    tNw, tNk = tFT2d.shape
    # find angle according to Kernel estimation
    try:
        grad, angle = getdopgrad(tFT2d,norm=(dk*vA/dw)) # normalised to 1/vA
        dopangle = (-const.PI/2-angle)# *180/const.PI # change direction in which angle is measured
        if dopangle > 0 or dopangle < maxang*const.PI/180:
            raise SystemError
    except:
        dopangle = defangle
    print('Kernel extracted angle : ',dopangle*180/const.PI)
    # setup plot
    fig,ax=plt.subplots(figsize=(10,8),nrows=3,layout='constrained')
    # plot FT2d
    ax[0].imshow(np.log10(tFT2d),origin='lower',interpolation='none',aspect='auto',cmap='magma',extent=extents)
    # ax[0].plot([0,100],[0,100],linestyle='--',color='k') # Alfven line
    ax[0].set_ylim(0,wmax_prime) ; ax[0].set_xlim(0,kmax_prime)
    ax[0].set_ylabel('Frequency '+r'$[\Omega_p]$',**tnrfont)
    ax[0].set_xlabel('Wavenumber '+r'$[\Omega_p/v_A]$',**tnrfont)
    # loop through shifted prime frequencies and find power
    freqs_prime = np.linspace(0,wmax_prime,tNw) # normalised units
    ynarr = freqs_prime*tNw/wmax_prime # np.linspace(0,tNw,tNw) # pixel units
    # plot harmonics
    for i in range(wmax_prime):
        ax[1].axvline(i,linestyle='--',color='darkgrey')
        ax[2].axvline(i,linestyle='--',color='darkgrey')
    # load & plot normal power spectra (dopangle = -90deg)
    omegas_power = read_pkl('omegas_power')/wcmin
    log10_psd = read_pkl('log10_power') - np.log10(dw) # psd=power/dw # log(psd)=log(power)-log(dw)
    ax[1].plot(omegas_power[omegas_power<wmax_prime],log10_psd[omegas_power<wmax_prime],color='k')
    ax[1].annotate(r'$\omega$',xy=(0.95,0.05),xycoords='axes fraction',va='bottom',ha='right',**tnrfont,color='k')

    # plot example lines
    ax[0].plot([0,kmax_prime],[5,5],linestyle='--',color='k')
    ax[0].plot([0,kmax_prime],[5,(1/np.tan(dopangle))*kmax_prime+5],linestyle='--',color='b')

    # find power along line
    power_line = getPowerLine(tFT2d,ynarr,freqs_prime,kmax_prime=kmax_prime,wmax_prime=wmax_prime,angle=dopangle)

    # formatting
    ax[0].annotate("{:.0f}".format(xi2),xy=(0.95,0.95),xycoords='axes fraction',**tnrfont,va='top',ha='right')
    ax[2].plot(freqs_prime,np.log10(power_line)-np.log10(dw),color='b')
    ax[2].annotate(r'$\omega^\prime$',xy=(0.95,0.05),xycoords='axes fraction',va='bottom',ha='right',**tnrfont,color='b')
    ax[2].set_xlabel('Frequency '+r'$[\Omega_p]$',**tnrfont)
    ax[1].set_ylabel(r'$\log_{10}($'+'PSD'+r'$)$',**tnrfont)
    ax[2].set_ylabel(r'$\log_{10}($'+'PSD'+r'$)$',**tnrfont)
    fig.savefig(home+'power_doppler_xi2_{:.0f}_ang_{:.2f}.png'.format(xi2,dopangle*180/const.PI),bbox_inches='tight')
    pass

def plotdopPower(home,sims,labels,minspec='Protons',wkmax=[20,45],xlims=[0,20],width=10,height=5,\
                omegalabel=True,freqlabel=False,leg=False):
    if freqlabel == omegalabel:
        omegalabel = True ; freqlabel = False

    # colors for each sim
    colors = plt.cm.rainbow(np.linspace(0,1,len(sims)))
    # freq plotting limits
    xmin,xmax=xlims
    fig,ax = plt.subplots(figsize=(width,height))
    for j in range(xmin,xmax):
        ax.axvline(j,color='darkgrey',linestyle='--')
    i=0
    for sim in sims:
        # load and setup
        simloc = getSimulation(sim)
        d0 = sdfread(0)
        # load FFT2d and plot
        FT2d = read_pkl('FT_2d_Magnetic_Field_Bz')
        Nw, Nk = FT2d.shape
        times = read_pkl('times')
        klim, wlim = getDispersionlimits(sim)
        kmin, kmax = klim
        wmin, wmax = wlim
        # frequency and Alfven speed
        wcmin = getCyclotronFreq(d0,minspec)
        vA = getAlfvenVel(d0)
        # normalised cell widths
        dk = (2*const.PI/getGridlen(d0))
        dw = (2*const.PI/times[-1])
        # extents = [0,kmax*vA/wcmin,0,wmax/wcmin]
        # thresh FT2d
        wmax_prime, kmax_prime = wkmax
        extents = [0,kmax_prime,0,wmax_prime]
        # threshold FT2d
        kthresh = int(FT2d.shape[1]*extents[1]*wcmin/vA/kmax)
        wthresh = int(FT2d.shape[0]*extents[-1]*wcmin/wmax)
        tFT2d = FT2d[:wthresh,:kthresh]
        tNw, tNk = tFT2d.shape

        # get theory angle
        Emin = 14.68
        dsth = -np.abs(li.getdoptheory(d0,Emin,wcmin,minspec,vA)) # make sure -ve
        dopangle = np.arctan(1/(dsth*dk*vA/dw)) # already correct direction
        print('dsth : ',dsth,'\ndopangle :',dopangle*180/const.PI)

        # # get kernel angle
        # grad, angle = li.getdopgrad(tFT2d,norm=(dk*vA/dw)) # grad normalised to 1/vA
        # dopangle = (-const.PI/2-angle) # change direction in which angle is measured
        # print('dgrad :',grad,'\ndopangle :',dopangle*180/const.PI)

        # plot doppler
        freqs_prime = np.linspace(0,wmax_prime,tNw) # normalised units
        thresh = (freqs_prime > xmin) & (freqs_prime < xmax)
        ynarr = freqs_prime*tNw/wmax_prime # np.linspace(0,tNw,tNw) # pixel units
        power_line = getPowerLine(tFT2d,ynarr,freqs_prime,kmax_prime=kmax_prime,wmax_prime=wmax_prime,angle=dopangle)
        ax.plot(freqs_prime,np.log10(power_line)-np.log10(dw),color=colors[i],label=labels[i])
    
        i+=1

    # zoom and formatting
    if omegalabel:
        ax.set_xlabel(r'$\omega^\prime/$'+getOmegaLabel(minspec),**tnrfont)	
    if freqlabel:
        ax.set_xlabel(r'$f^\prime/$'+getFreqLabel(minspec),**tnrfont)	
    if leg:
        ax.legend(loc='best',labelspacing=0.1,borderpad=0.1,ncol=1)
    ax.set_xlim(xmin,xmax)
    ax.set_ylabel('PSD',**tnrfont)
    ax.locator_params(axis='x',nbins=6)
    # plt.show()
    fig.savefig(home+'power_compare_dopw_{}_{}.png'.format(xmin,xmax),bbox_inches='tight')
    
    return None
    

if __name__=='__main__':
    home = '/storage/space2/phrmsf/lowres_D_He3/'
    # get sims
    sims = np.sort([i for i in os.listdir(home) if 'p_90' in i])
    xiHe3 = np.array([int(i[2:4]) for i in sims])
    sims = np.array([home+i for i in sims])
    # sims = [home+'0_00_p_90']
    # xiHe3= [0]
    plotdopPower(home,sims[1:],xiHe3[1:],xlims=[12,20],leg=False,height=3)
    sys.exit()
    for i in range(len(sims)):
        print(sims[i])
        plotPowerLine(home,sims[i],xiHe3[i])
    
    # find vdop via kernel estimation (+ plot)

    # extract values along line, normalise to No. cells in line

    # plot power spectra equivalent along with standard power spectra



