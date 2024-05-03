
from func_load import *
import scipy


def getdopgrad(FT2d,logthresh=1.8,norm=1,kernel='custom',dydxrange=(-0.5,0.5)):
    """
        Get doppler gradient in range of dydx space
        Get this first from angles, then transfer to angle in terms of 
        real coords as angles are noisy
    """
    # make copy of FT2d and thresh
    tFT2d = FT2d.copy()
    tFT2d[FT2d < logthresh] = 0
    tFT2d[FT2d > logthresh] = 1
    
    # Kernel gradient map
    _,kGangle = Kernel(tFT2d,kernel=kernel) # list_new func : kernel = 'scharr' or 'sobel'
    # plt.imshow(kGangle,**imkwargs) ; plt.show()

    # remove edge abberations
    kGangle = kGangle[1:-1,1:-1]

    # flatten to take hist
    kGangle = kGangle.flatten()

    # convert to all negative angles (easier to calc real gradient)
    kGangle[kGangle>0] -= const.PI # phi = pi - theta
    
    # # remove gradients larger than 90deg (or limiting angle)
    # kGangle = kGangle[np.abs(kGangle)<const.PI/2]

    # convert grad from angle to velocity (convert +-inf to nan)
    dydx = np.nan_to_num(np.tan(kGangle),posinf=np.nan,neginf=np.nan)
    
    # remove nan & zero values
    dydx = dydx[~np.isnan(dydx)]
    dydx = dydx[dydx!=0]

    # normalise and get hist
    dydx/=norm
    counts,bins=np.histogram(dydx,bins=1000,range=dydxrange,density=True)
    # plt.hist(dydx,bins=1000,range=dydxrange,density=True) ; plt.show()
    
    # only get negative gradients
    # get grad and calc corresponding angle
    grad = -np.abs(bins[np.argmax(counts)]) # negative gradient
    angle = np.arctan(grad*norm)
    print(grad,angle*180/const.PI)

    return grad, angle # radians 

def getPowerLine(FT2d,ynarr,freqs,kmax_prime=100,wmax_prime=20,angle=-90):
    Nw, Nk = FT2d.shape
    # find power along line
    power_line = np.zeros(len(freqs))
    lsizearr = []
    for i in range(len(ynarr)):
        yn = ynarr[i]
        xlim,ylim = ds.LineBoxIntersection(Ys=Nw,Ye=0,Xn=0,Yn=yn,XL=Nk,YL=Nw,\
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

def plotPowerLine(home,sim,xi2,kmax_prime=20,wmax_prime=10,minspec='Protons',angle=-85):
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
    extents = [0,kmax_prime,0,wmax_prime]
    # threshold FT2d
    kthresh = int(FT2d.shape[1]*extents[1]*wcmin/vA/kmax)
    wthresh = int(FT2d.shape[0]*extents[-1]*wcmin/wmax)
    tFT2d = FT2d[:wthresh,:kthresh]
    tNw, tNk = tFT2d.shape
    # find angle according to Kernel estimation
    grad, angle = getdopgrad(tFT2d,norm=(dk*vA/dw)) # normalised to 1/vA
    dopangle = (-const.PI/2-angle)*180/const.PI
    print(dopangle)
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
    log10_psd = read_pkl('log10_power') - np.log10(dw)
    ax[1].plot(omegas_power[omegas_power<wmax_prime],log10_psd[omegas_power<wmax_prime],color='k')
    ax[1].annotate(r'$\omega$',xy=(0.95,0.05),xycoords='axes fraction',va='bottom',ha='right',**tnrfont,color='k')

    # plot example lines
    ax[0].plot([0,kmax_prime],[5,5],linestyle='--',color='k')
    ax[0].plot([0,kmax_prime],[5,(1/np.tan(dopangle*const.PI/180))*kmax_prime+5],linestyle='--',color='b')

    # find power along line
    power_line = getPowerLine(tFT2d,ynarr,freqs_prime,kmax_prime=kmax_prime,wmax_prime=wmax_prime,angle=dopangle)

    # formatting
    ax[0].annotate("{:.0f}".format(xi2),xy=(0.95,0.95),xycoords='axes fraction',**tnrfont,va='top',ha='right')
    ax[2].plot(freqs_prime,np.log10(power_line)-np.log10(dw),color='b')
    ax[2].annotate(r'$\omega^\prime$',xy=(0.95,0.05),xycoords='axes fraction',va='bottom',ha='right',**tnrfont,color='b')
    ax[2].set_xlabel('Frequency '+r'$[\Omega_p]$',**tnrfont)
    ax[1].set_ylabel(r'$\log_{10}($'+'PSD'+r'$)$',**tnrfont)
    ax[2].set_ylabel(r'$\log_{10}($'+'PSD'+r'$)$',**tnrfont)
    fig.savefig(home+'power_doppler_xi2_{:.0f}_ang_{:.2f}.png'.format(xi2,dopangle),bbox_inches='tight')
    pass

if __name__=='__main__':
    import dopplerShift.line_integrate as ds
    home = '/storage/space2/phrmsf/lowres_D_He3/'
    # get sims
    sims = np.sort([i for i in os.listdir(home) if 'p_90' in i])
    xiHe3 = np.array([int(i[2:4]) for i in sims])
    sims = np.array([home+i for i in sims])
    # sims = [home+'0_00_p_90']
    # xiHe3= [0]
    for i in range(len(sims)):
        print(sims[i])
        plotPowerLine(home,sims[i],xiHe3[i])
    
    # find vdop via kernel estimation (+ plot)

    # extract values along line, normalise to No. cells in line

    # plot power spectra equivalent along with standard power spectra



