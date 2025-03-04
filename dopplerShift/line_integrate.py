
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

def getdoptheory(file0,wcmin,minspec,vA,Emin=14.68,uperp_vA=0.9):
    """
    DESCRIPTION HERE
        IN:
            file0       : sdf read object, 0th file (00000.sdf)
            wcmin       : float, cyclotron frequency of minority species
            minspec     : str, name of minority species
            vA          : float, Alfven speed [m/s]
            Emin        : float, energy of minority species [MeV]
            uperp_vA    : float, the ratio between the perpendicular birth velocity of minspec and the Alfven wave speed vA
        OUT:
            dsth        : float, returns the Doppler velocity as predicted by theory using pitch-angle and magnetic angle
    """
    # get minority birth uperp (in terms of vA)
    try:
        theta_B,_ = getMagneticAngle(file0)
    except:
        theta_B = 89*np.pi/180
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
        IN:
            Ys, Ye  : int, start and end point of line on y-axis
            Yn, Xn  : int, offset of line mid-point
            YL, XL  : int, total length of box in Y and X (pixel coordinates)
            theta   : float, positive angle clockwise from north [rad], represents 1/gradient of line
        OUT:
            xlim    : arr, array of limits of box in x direction
            ylim    : arr, array of limits of box in y direction
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

def getPowerLine(FT2d,ynarr,freqs,kmax_prime=100,wmax_prime=20,angle=-3.14/2):
    """
    DESCRIPTION HERE
        IN:
            param1          :
        OUT:

    """
    Nw, Nk = FT2d.shape
    # find power along line
    power_line = np.zeros(len(freqs))
    lsizearr = []
    for i in range(len(ynarr)):
        yn = ynarr[i]
        xlim,ylim = LineBoxIntersection(Ys=Nw,Ye=0,Xn=0,Yn=yn,XL=Nk,YL=Nw,\
                                            theta=angle)
        # number of cells within line
        lsizearr.append(np.around(np.sqrt(((xlim[1]-xlim[0]))**2 + ((ylim[1]-ylim[0]))**2))) # No. pixels
        # pixel coordinate points along line
        x = np.linspace(xlim[0],xlim[1],int(lsizearr[-1]))
        y = np.linspace(ylim[0],ylim[1],int(lsizearr[-1]))
        # values along line
        zi = scipy.ndimage.map_coordinates(FT2d,np.vstack((y,x)),order=0)
        # power along line
        power_line[i] = np.sum(zi**2) #/len(zi)
    return power_line
