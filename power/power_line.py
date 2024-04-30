
from func_load import *
import scipy

def getdopangle(FT2d,logthresh=1.8):
    tFT2d = FT2d.copy()
    tFT2d[FT2d < 1.8] = 0
    tFT2d[FT2d > 1.8] = 1
    # Kernel gradient map
    _,kGangle = Kernel(tFT2d,kernel='sobel') # list_new func : kernel = 'scharr' or 'sobel'
    kGangle = kGangle.flatten()
    # convert to all negative angles (easier to calc real gradient)
    kGangle[kGangle>0] -= const.PI

    # convert grad from angle to velocity (convert +-inf to nan)
    dydx = np.nan_to_num(np.tan(kGangle),posinf=np.nan,neginf=np.nan)
    # remove nan & zero values
    dydx = dydx[~np.isnan(dydx)]
    dydx = dydx[dydx!=0]
    return dydx

if __name__=='__main__':
    import dopplerShift.line_integrate as ds
    home = '/storage/space2/phrmsf/lowres_D_He3/'
    # get sims
    sims = [home+'0_00_p_90']
    for sim in sims:
        simloc = getSimulation(sim)
        # load FFT2d and plot
        FT2d = read_pkl('FT_2d_Magnetic_Field_Bz')
        Nw, Nk = FT2d.shape
        times = read_pkl('times')
        klim, wlim = getDispersionlimits(sim)
        kmin, kmax = klim
        wmin, wmax = wlim
        print(kmin,kmax,wmin,wmax)
        wcmin = getCyclotronFreq(sdfread(0),'Protons')
        vA = getAlfvenVel(sdfread(0))
        # normalised cell widths
        dk = (kmax/Nk)*vA/wcmin
        dw = (wmax/Nw)/wcmin
        # extents = [0,kmax*vA/wcmin,0,wmax/wcmin]
        # thresh FT2d
        kmax_prime = 50 # 80 # kmax*vA/wcmin
        wmax_prime = 20
        extents = [0,kmax_prime,0,wmax_prime]
        kthresh = int(FT2d.shape[1]*extents[1]*wcmin/vA/kmax)
        wthresh = int(FT2d.shape[0]*extents[-1]*wcmin/wmax)
        tFT2d = FT2d[:wthresh,:kthresh]
        tNw, tNk = tFT2d.shape
        print(tFT2d.shape)
        # setup plot
        fig,ax=plt.subplots(figsize=(10,8),nrows=3,layout='constrained')
        # plot FT2d
        ax[0].imshow(np.log10(tFT2d),origin='lower',interpolation='none',aspect='auto',cmap='magma',extent=extents)
        ax[0].plot([0,100],[0,100],linestyle='--',color='k')
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
        log10_power = read_pkl('log10_power')
        ax[1].plot(omegas_power[omegas_power<wmax_prime],log10_power[omegas_power<wmax_prime],color='b')
        ax[1].annotate(r'$\omega$',xy=(0.95,0.05),xycoords='axes fraction',va='bottom',ha='right',**tnrfont,color='b')
        # find angle according to Kernel estimation
        # threshold FT2d
        dydx = getdopangle(tFT2d)*(dw/dk)
        fig.clf()
        plt.hist(dydx,bins=100) ; plt.show()
        # plot example lines
        ax[0].plot([0,kmax_prime],[5,5],linestyle='--',color='b')
        ax[0].plot([0,kmax_prime],[5,(1/np.tan(dopangle*const.PI/180))*kmax_prime+5],linestyle='--',color='r')
        # find power along line
        power_line = np.zeros(len(freqs_prime))
        lsizearr = []
        for i in range(len(ynarr)):
            yn = ynarr[i]
            xlim,ylim = ds.LineBoxIntersection(Ys=tNw,Ye=0,Xn=0,Yn=yn,XL=tNk,YL=tNw,\
                                                theta=dopangle)
            # ax[0].plot(kmax_prime*xlim/tNk,wmax_prime*ylim/tNw,color='k') # pixel coordinates
            # number of cells within line
            lsizearr.append(np.around(np.sqrt(((xlim[1]-xlim[0]))**2 + ((ylim[1]-ylim[0]))**2))) # No. pixels
            # pixel coordinate points along line
            x = np.linspace(xlim[0],xlim[1],int(lsizearr[-1]))
            y = np.linspace(ylim[0],ylim[1],int(lsizearr[-1]))
            # values along line
            zi = scipy.ndimage.map_coordinates(tFT2d,np.vstack((y,x)),order=0)
            # power along line
            power_line[i] = np.sum(zi**2)
        ax[2].plot(freqs_prime,np.log10(power_line),color='r')
        ax[2].annotate(r'$\omega^\prime$',xy=(0.95,0.05),xycoords='axes fraction',va='bottom',ha='right',**tnrfont,color='r')
        ax[1].set_ylim(5)
        ax[2].set_ylim(5)
        ax[2].set_xlabel('Frequency '+r'$[\Omega_p]$',**tnrfont)
        ax[1].set_ylabel(r'$\log_{10}($'+'Power'+r'$)$',**tnrfont)
        ax[2].set_ylabel(r'$\log_{10}($'+'Power'+r'$)$',**tnrfont)
        plt.show()

    # find vdop via kernel estimation (+ plot)

    # extract values along line, normalise to No. cells in line

    # plot power spectra equivalent along with standard power spectra



