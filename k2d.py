
def make2D(rowval,colval,val,rowlim=(None,None),collim=(None,None),bins=(1000,1000),limits=False):
    if rowval.shape != colval.shape and rowval.shape != val.shape: # make sure same shape
            raise SystemError
    if rowlim[0] != None or collim[0] != None: # check if limits applied
        # thresh to limit size
        thresh = (rowlim[0]<rowval) & (rowval<rowlim[1]) & (collim[0]<colval) & (colval<collim[1]) 
        rowval = rowval[thresh] ; colval = colval[thresh] ; val = val[thresh]
        # min, max values from data
        rowmin, rowmax = [np.min(rowval),np.max(rowval)]
        colmin, colmax = [np.min(colval),np.max(colval)]
    if bins[0] == None: # no bins
        # unique values
        urow,urowind = np.unique(rowval,return_index=True)
        ucol,ucolind = np.unique(colval,return_index=True)
        nx,ny=[len(urow),len(ucol)]
    else: # set bin size
        nx,ny=bins
        # min max values from user input
        rowmin, rowmax = rowlim
        colmin, colmax = collim

    # arrays between max and min values
    rowarr = np.linspace(rowmin,rowmax,nx)
    colarr = np.linspace(colmin,colmax,ny)
    Z = np.zeros((len(rowarr),len(colarr)))
    for k in range(len(rowval)):
        print(k)
        i = np.where(rowval[k] >= rowarr)[0][-1] # last index corresponding to row
        j = np.where(colval[k] >= colarr)[0][-1] # last index corresponding to column
        if Z[i,j] < val[k]:
            Z[i,j]=val[k] # assign highest growth rate

    if limits:
        return Z, [colmin,colmax,rowmin,rowmax]
    else:
        return Z



if __name__=='__main__':
    from func_load import *
    import lib.frequencies as freq
    import matplotlib.tri as tri
    B0 = 2.07
    n0 = 1.7e19
    xialpha = 1.5e-4
    xiT = 0.11
    xiD = (1)*(1-(1)*xiT-2*xialpha)
    # get all frequencies
    wce = freq.CyclotronFreq(B0,1,const.me)
    wcD = freq.CyclotronFreq(B0,1,getMass('Deuterons'))
    wcT = freq.CyclotronFreq(B0,1,getMass('Tritons'))
    wcalpha = freq.CyclotronFreq(B0,2,getMass('Alphas'))
    wpe = freq.PlasmaFreq(n0,1,const.me)
    wpD = freq.PlasmaFreq(n0*xiD,1,getMass('Deuterons'))
    wpT = freq.PlasmaFreq(n0*xiT,1,getMass('Tritons'))
    wpalpha = freq.PlasmaFreq(n0*xialpha,2,getMass('Alphas'))
    # normalisation k
    mD = getMass('Deuterons')
    mT = getMass('Tritons')
    malpha = getMass('Alphas')
    vA = B0/np.sqrt(const.mu0*n0*(xiD*mD+xiT*mT+xialpha*malpha))
    k0 = wcalpha/vA

    # frequencies and arrays
    dw = 0.1
    omegas = wcalpha * np.arange(0,15,dw)
    wpf = [wpe,wpD,wpT,wpalpha]
    wcf = [wce,wcD,wcT,wcalpha]
    dk = 0.025
    dtheta = dk#/(omegas[-1]/wcalpha)
    theta = np.arange(0,90.1,dtheta)*const.PI/180
    # theta = [89*const.PI/180]

    # assign empty arrays    
    kpara = np.zeros((len(theta)*len(omegas)))
    kperp = np.zeros((len(theta)*len(omegas)))
    omegas_arr = np.zeros((len(theta)*len(omegas)))

    kperpmax = 20
    kparamax = 4
    bins = 200
    try:
        # raise SystemError
        Z = read_pkl('k2d_kperp_{}_kpara_{}'.format(kperpmax,kparamax))
    except:
        # loop over angles
        for i in range(len(theta)):
            _,k2,_ = coldplasmadispersion_analytical(omegas,wpf=wpf,wcf=wcf,theta=theta[i])
            # assign
            kpara[i*len(omegas):(i+1)*len(omegas)] = k2*np.cos(theta[i])
            kperp[i*len(omegas):(i+1)*len(omegas)] = k2*np.sin(theta[i])
            omegas_arr[i*len(omegas):(i+1)*len(omegas)] = omegas
            # plt.scatter(k2*np.sin(theta[i]),k2*np.cos(theta[i]),c=omegas,edgecolor='none')
            # plt.annotate(theta[i],(k2[-1]*np.sin(theta[i]),k2[-1]*np.cos(theta[i])))

        # plt.scatter(kperp,kpara,c=omegas_arr/wcalpha,edgecolor='none',cmap='Accent')
        # plt.show()
        print(kpara.shape,kperp.shape,omegas_arr.shape)
        Z = make2D(kpara,kperp,omegas_arr/wcalpha,rowlim=(0,kparamax*k0),collim=(0,kperpmax*k0),limits=False,bins=(bins,bins))
        dumpfiles(Z,'k2d_kperp_{}_kpara_{}'.format(kperpmax,kparamax))

    fig,ax=plt.subplots(nrows=3,figsize=(8,15),sharex=True)
    fig.subplots_adjust(hspace=0.075)
    im = ax[0].imshow(Z,origin='lower',aspect='auto',cmap='jet',interpolation='none',extent=[0,kperpmax,0,kparamax],clim=(0,15))
    # cbar = plt.colorbar(im)
    cntr = ax[0].contour(Z,levels=np.arange(0,kperpmax+1,1),origin='lower',aspect='auto',colors='k',extent=[0,kperpmax,0,kparamax])
    plt.clabel(cntr, inline=1, fontsize=12)

    # compare to high freq equivalent
    tkperp = np.linspace(0,kperpmax,bins)*k0
    tkpara = np.linspace(0,kparamax,bins)*k0
    KPERP,KPARA=np.meshgrid(tkperp,tkpara)
    K2 = KPARA**2 + KPERP**2
    W2 = ((vA**2)/2)*(K2 + KPARA**2 + (K2*KPARA**2)*((vA**2)/(wcalpha**2)) + \
            ((K2 + KPARA**2 + (K2*KPARA**2)*(vA**2)/(wcalpha**2))**2 - 4*K2*KPARA**2)**0.5)
    ax[1].imshow(np.sqrt(W2)/wcalpha,origin='lower',aspect='auto',cmap='jet',interpolation='none',extent=[0,kperpmax,0,kparamax],clim=(0,15))
    cntr = ax[1].contour(np.sqrt(W2)/wcalpha,levels=np.arange(0,kperpmax+1,1),origin='lower',colors='w',extent=[0,kperpmax,0,kparamax])
    plt.clabel(cntr, inline=1, fontsize=12)
    # ratio between solutions 
    ax[2].imshow((Z*wcalpha)/np.sqrt(W2),origin='lower',aspect='auto',cmap='jet',interpolation='none',extent=[0,kperpmax,0,kparamax])
    cntr = ax[2].contour((Z*wcalpha)/np.sqrt(W2),levels=np.linspace(1,2,5),origin='lower',colors='w',extent=[0,kperpmax,0,kparamax])    
    plt.clabel(cntr, inline=1, fontsize=12)
    # labels and annotations
    fig.supylabel('Parallel Wavenumber' + '  '+r'$[\Omega_i/V_A]$',**tnrfont)
    fig.supxlabel('Perpendicular Wavenumber' + '  '+r'$[\Omega_i/V_A]$',**tnrfont)
    ax[0].annotate('FAW Stix',xy=(0.9,0.9),xycoords='axes fraction',ha='right',va='bottom',color='w',**tnrfont)
    ax[1].annotate('FAW Cook',xy=(0.9,0.9),xycoords='axes fraction',ha='right',va='bottom',color='w',**tnrfont)
    ax[2].annotate('Stix/Cook',xy=(0.9,0.9),xycoords='axes fraction',ha='right',va='bottom',color='w',**tnrfont)
    # colorbar
    p0 = ax[0].get_position().get_points().flatten()
    p1 = ax[-1].get_position().get_points().flatten()
    cbar = fig.add_axes([p0[2]+0.02, p1[1], 0.01, p0[3]-p1[1]]) # [left bottom width height]
    clbar = plt.colorbar(im, cax=cbar, orientation='vertical', spacing='proportional')
    clbar.ax.set_ylabel('$\omega/\Omega_i$',**tnrfont)
    plt.savefig('top_FAWStix_bottom_FAWJWSCook.png',bbox_inches='tight')
    plt.show()
    
