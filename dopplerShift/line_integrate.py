
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
    # print(grad,angle*180/const.PI)

    return grad, angle # radians 

def getdoptheory(file0,Emin,wcmin,minspec,vA,uperp_vA=0.9):
    """
        file0 ; sdf read object
        Emin ; MeV
        wcmin ; 1/s
        minspec ; name
        vA ; m/s
    """
    # get minority birth uperp (in terms of vA)
    theta_B,_ = getMagneticAngle(file0)
    # get minority birth energy, mass and velocity
    umin = np.sqrt(2*Emin*1e6*const.qe/getMass(minspec))
    uperp = uperp_vA*vA
    # get pitch angle
    pitch_angle = np.arcsin(uperp/umin)
    # get dsth
    dsth = (umin/vA)*np.cos(theta_B)*np.cos(pitch_angle)
    return dsth # normalised to vA

def LineBoxIntersection(Ys,Ye,Xn,Yn,XL,YL,theta):
    """
    Find intersection points between a line with eqn (Y-Yn)=(X-Xn)/tan(theta)
        In:
            Ys, Ye : start and end point of line on y-axis
            Yn, Xn : offset of line mid-point
            YL, XL : total length of box in Y and X (pixel coordinates)
            theta : positive angle clockwise from north [rad], represents 1/gradient of line
    """
    # +ve angle clockwise from north
    Xs = (Ys-Yn)*np.tan(theta)+Xn ; Xe = (Ye-Yn)*np.tan(theta)+Xn
    # check initial co-ordinates, if wrong then use other bounds 
    if Xs < 0:
        Xs = 0
        Ys = Yn+((Xs-Xn)/np.tan(theta))
    if Xs > XL:
        Xs = XL
        Ys = Yn+((Xs-Xn)/np.tan(theta))
    if Xe > XL:
        Xe = XL
        Ye = Yn+((Xe-Xn)/np.tan(theta)) 
    if Xe < 0:
        Xe = 0
        Ye = Yn+((Xe-Xn)/np.tan(theta)) 
    xlim = [Xs,Xe] ; ylim = [Ys,Ye]
    return np.array(xlim), np.array(ylim)

def line_integrate(Z,xposini=[],yposini=[],angles=[],rowlim=(-4,4),collim=(0,15),norm=[1,1],lsize=None,label='',colorarr=True):
    """
        In:
            Z : 2d matrix of data points to calculate integral over
            x/yposini : x/y positions, in units of x/y-axis px coords, as to where centre line integral
            angles : angles over which to loop, defaults to 0 -> 90deg
        Out:
            intensity : summated line intensity normalised to the number of cells within the line (sum(zi)/len(zi))
            dopangmax : array of angles per yposini that correspond to the strongest total intensity along the line
            Zij       : 3d array of the line intensity [shape of (len(yposini), len(angles), len(zi))]
    """
    Ny, Nx = Z.shape
    if not lsize:
        lsize = np.sqrt(Nx**2+Ny**2) # max length of line within imshow square (px coords)
    # xposini and yposini same len
    if xposini.shape != yposini.shape:
        xposini = [0] ; yposini = [0]
    # set angle array if empty
    if angles==[]:
        angles = np.linspace(0,-180,360)
    # color array
    if colorarr:
        colors = plt.cm.rainbow(np.linspace(0,1,len(angles)))
    else:
        colors = ['k']*len(angles)

    # setup empty arrays
    dopmaxang=np.zeros(len(yposini))
    zi = np.zeros((len(angles),int(lsize)))
    intensity = np.zeros((len(yposini),len(angles)))
    rangles = np.zeros(len(angles)) # real angles
    Zij = np.zeros((len(yposini),len(angles),int(lsize)))
    
    # limits of box, normalised units
    rowlim = np.array(rowlim) ; collim = np.array(collim)
    xmin, xmax = (collim/norm[0])
    ymin, ymax = (rowlim/norm[1])
    extents = [xmin,xmax,ymin,ymax]
    
    # transform to pixel coords
    px_xposini = (np.array(xposini)/(xmax-xmin))*Nx
    py_yposini = (np.array(yposini)/(ymax-ymin))*Ny

    # loop through angles
    for i in range(len(py_yposini)):
        Xn = px_xposini[i] # number of cells in x
        Yn = py_yposini[i] # number of cells in y
        # Xn = px_xposini[i]
        # Yn = Ny/2 # centred at centre of image (useful for symmetric plots)
        for j in range(len(angles)):
            # find intersection points between line and box (pixel coords)
            Ystart = Ny ; Yend = 0
            xlim, ylim = LineBoxIntersection(Ystart,Yend,Xn,Yn,Nx,Ny,angles[j]*180/np.pi)
            # lsize = np.sqrt((xlim[-1]-xlim[0])**2+(ylim[-1]-ylim[0])**2)
            # find data points along line
            x = np.linspace(xlim[0],xlim[1],int(lsize))
            y = np.linspace(ylim[0],ylim[1],int(lsize))
            zi[j,:] = scipy.ndimage.map_coordinates(Z/norm[0],np.vstack((y,x))) # normalised growth rates
            Zij[i,j,:]=zi[j,:]
            # convert all to real coordinates
            print(xlim,ylim)
            xlim = (collim[0]/norm[0])+xlim*((collim[1]-collim[0])/norm[0])/Nx
            ylim = (rowlim[0]/norm[1])+ylim*((rowlim[1]-rowlim[0])/norm[1])/Ny
            rangles[j] = np.arctan((xlim[-1]-xlim[0])/(ylim[-1]-ylim[0]))*180/np.pi # [deg]
            # summate all points # "growth per (dkpara,dw) cell" 
            tzi = zi[j,:]
            tzi[tzi < 0] = 0 # no negative growth rates affecting intensity extraction
            intensity[i,j] = np.sum(tzi**2)#/len(tzi) # normalise to number of cells along line
            # # example plot
            # if plot:
            #     ax[0].imshow(Z,origin='lower',interpolation='none',aspect='auto',cmap='summer',extent=extents)
            #     ax[0].plot(xlim,ylim,color=colors[j],linestyle='--')#,alpha=1/len(angles))
            #     ax[1].plot(np.linspace(xlim[0],xlim[1],len(zi[j,:])),zi[j,:],color=colors[j],label=rangles[j])#,alpha=1/len(angles))

        # find maximum intensity as a function of angle per yposini
        maxintarg = np.argmax(intensity[i,:])
        maxang = rangles[maxintarg]
        print('Max intensity angle [deg]: ',maxang)
        dopmaxang[i] = maxang
        # # example plot
        # ax[0].set_ylim(tparamin,yparamax) ; ax[0].set_xlim(xmin,xmax)
        # # ax[0].set_xlabel('Frequency '+r'$[\Omega_i]$',**tnrfont)
        # ax[0].set_ylabel('Parallel Wavenumber '+r'$[\Omega_i/V_A]$',**tnrfont)
        # ax[1].set_xlabel('Frequency '+r'$[\Omega_i]$',**tnrfont)
        # ax[1].set_ylabel('Growth Rate '+r'$[\Omega_i]$',**tnrfont)
        # ax[1].legend(loc='best')
        # fig.savefig('all_lines_XI2_{}.png'.format(label))
    return intensity, dopmaxang, Zij