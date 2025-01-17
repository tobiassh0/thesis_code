
from func_load import *
import dopplerShift.line_integrate as li


def getIntegrateDoppler(sims,FT2darr,normspecies='Protons',wkmax=[20,20],logthresh=1.8):
	"""
	Extracts the maximum point in a smaller array (so don't use edge of array and bias sample) of the 2d FFT
	the loops through multiple angles (angle +ve clockwise from North) and extracts the values of the threshed
	2d FFT array, summates them then plots this integrand vs its corresponding angle. Can then find the angle
	of doppler shift and convert this to a velocity, plotting this atop the 2d FFT for multiple l 
	harmonics.  
		IN:
			sims:			list of sims to analyse and extract each gradient from
			FT2darr:		list of 2d FFTs corresponding to the sims in list sims 
			normspecies:the normalisation species used for w,k space
			wkmax:		the limits of the 2d FFT [wmax,kmax] which to plot, in units of normspecies normalisation 
			logthresh:	the log value of the thresh which to apply to the 2d FFT 
		OUT:
			Plots the gradient angle vs. normalised integral to the number of cells np.sum(zi)/len(zi) 
	"""	
	for i in range(len(sims)):
		# setup
		sim_loc = getSimulation(sim)	
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
		swmax = (wkmax[0]-10)*wnorm ; skmax = (wkmax[1]-10)*knorm
		# reduce limit and find max
		sFT2d = np.log10(FT2d[:int(nw*swmax/wlim),:int(nk*skmax/klim)])
		FT2d = np.log10(FT2d[:int(nw*wmax/wlim),:int(nk*kmax/klim)])
		(nw,nk) = FT2d.shape

		# threshold FT2d
		tFT2d = FT2d.copy() ; tsFT2d = sFT2d.copy() # copy arr
		tFT2d[FT2d < logthresh] = 0
		tFT2d[FT2d > logthresh] = 1
		tsFT2d[sFT2d < logthresh] = 0
		tsFT2d[sFT2d > logthresh] = 1
		FT2d = tFT2d # replace old FT2d 
		sFT2d = tsFT2d # replace old FT2d 
		del tFT2d, tsFT2d # destroy temp arr

		# extract maximum point along FAW
		argind = np.unravel_index(sFT2d.argmax(),sFT2d.shape)
		ym,xm = argind[0]*swmax/sFT2d.shape[0], argind[1]*skmax/sFT2d.shape[1]
		del sFT2d # remove smaller matrix
		print(xm/knorm,ym/wnorm)
		angles = np.linspace(const.PI/2,const.PI,1000) # between -180deg and -90deg
		# varying angle of line
		integ = np.zeros(len(angles))
		for i in range(len(angles)):
			# find image edge intercepts
			xi = (xm/knorm-(ym/wnorm)/np.tan(angles[i])) # algebra
			yi = (ym/wnorm-(xm/knorm)*np.tan(angles[i])) # normalised
			# ax.plot([0,xi],[yi,0],linestyle='--',color='white')
			# convert to pixel coordinates
			xi*=nk/(wkmax[1]); yi*=nw/(wkmax[0])
			x,y=np.linspace(0,xi,1000),np.linspace(0,yi,1000)
			# map coordinates to integrate
			zi = scipy.ndimage.map_coordinates(FT2d, np.vstack((y,x)))
			# integrate and append to array
			integ[i]=np.sum(zi/len(zi)) # normalise to integ per cell
		# maxangle (shared area)
		maxangle = angles[np.argmax(integ)]
		# plot angles vs integral
		fig,ax=plt.subplots(nrows=2,figsize=(6,8))
		ax[0].plot(angles*180/const.PI,integ)
		ax[0].plot([0,0],[0,maxangle],color='k',linestyle='--')
		ax[1].imshow(FT2d,**imkwargs,extent=[0,wkmax[1],0,wkmax[0]])
		# doppler shifted lines
		dsv = np.tan(maxangle) * (dw/dk)/vA		
		kx = np.linspace(0,wkmax[1],100)*knorm
		for l in range(0,int(wkmax[0]),1):
			w = wnorm*np.ones(len(kx))*l
			ww = w + (dsv*vA)*kx
			ax[1].plot(kx/knorm,ww/wnorm,color='white',linestyle='--')
		ax[1].set_xlabel(r'$kv_A/$'+getOmegaLabel(normspecies),**tnrfont)
		ax[1].set_ylabel(r'$\omega/$'+getOmegaLabel(normspecies),**tnrfont)
		ax[1].set_ylim(0,wkmax[0]) ; ax[1].set_xlim(0,wkmax[1])

		plt.show()
		sys.exit()
	return None	

def plotPowerLine(home,sim,xi2,wkmax=[10,20],minspec='Protons',maxang=-90,defangle=-85):
    """
    DESCRIPTION HERE
        IN:
            param1          :
        OUT:

    """
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
    power_line = li.getPowerLine(tFT2d,ynarr,freqs_prime,kmax_prime=kmax_prime,wmax_prime=wmax_prime,angle=dopangle)

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
    """
    DESCRIPTION HERE
        IN:
            param1          :
        OUT:

    """
    if freqlabel == omegalabel:
        omegalabel = True ; freqlabel = False

    # colors for each sim
    colors = plt.cm.rainbow(np.linspace(0,1,len(sims)))
    # freq plotting limits
    xmin,xmax=xlims
    fig,ax = plt.subplots(figsize=(width,height))
    for j in range(xmin,xmax):
        ax.axvline(j,color='darkgrey',linestyle='--')
    for i in range(len(sims)):
        # load and setup
        simloc = getSimulation(sims[i])
        d0 = sdfread(1) # 0
        # load FFT2d and plot
        FT2d = read_pkl('FT_2d_Magnetic_Field_Bz')
        Nw, Nk = FT2d.shape
        times = read_pkl('times')
        klim, wlim = getDispersionlimits(sims[i])
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

        # get doppler line power values
        freqs_prime = np.linspace(0,wmax_prime,tNw) # normalised units
        thresh = (freqs_prime > xmin) & (freqs_prime < xmax)
        ynarr = freqs_prime*tNw/wmax_prime # np.linspace(0,tNw,tNw) # pixel units
        power_line = li.getPowerLine(tFT2d,ynarr,freqs_prime,kmax_prime=kmax_prime,wmax_prime=wmax_prime,angle=dopangle)

        # plot PSD equivalent line
        ax.plot(freqs_prime,np.log10(power_line)-np.log10(dw),color=colors[i],label=labels[i])
    
    # frequency label (omega or f)
    if omegalabel:
        ax.set_xlabel(r'$\omega^\prime/$'+getOmegaLabel(minspec),**tnrfont)	
    if freqlabel:
        ax.set_xlabel(r'$f^\prime/$'+getFreqLabel(minspec),**tnrfont)	

    # legend
    if leg:
        ax.legend(loc='upper left',labelspacing=0.1,borderpad=0.1,columnspacing=0.1,ncol=len(sims))

    # zoom and formatting
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(-2,5)
    ax.set_ylabel('PSD',**tnrfont)
    ax.locator_params(axis='x',nbins=6)

    # plt.show()
    fig.savefig(home+'power_compare_dopw_{}_{}_short.png'.format(xmin,xmax),bbox_inches='tight')
    
    return None
    

if __name__=='__main__':
    # # D-He3
    # home = '/storage/space2/phrmsf/lowres_D_He3/'
    home = '/run/media/phrmsf/My Passport/simulations/D-He3/pklfiles-lowres_D_He3/'
    # get sims
    sims = np.sort([i for i in os.listdir(home) if 'p_90' in i])
    sims = sims[1:] # remove 0%
    sims = sims[::2] # every other simulation
    xiHe3 = np.array([int(i[2:4]) for i in sims])
    sims = np.array([home+i for i in sims])
    
    plotdopPower(home,sims,xiHe3,xlims=[12,20],leg=True,height=3)
    sys.exit()