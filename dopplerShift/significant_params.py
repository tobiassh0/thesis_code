
from func_load import *
from scipy import ndimage

def plot_significant_dop(l=4,Nx=1000,collim=(0,10),rowlim=(-1,1),Ny=1000,_levels=True,_Deltaspacing=0.1):

    extent=[collim[0],collim[1],rowlim[0],rowlim[1]]
    # setup grid
    uperp_vA = np.linspace(collim[0],collim[1],Nx)
    magnetic_field_angle = np.linspace(rowlim[0],rowlim[1],Ny)
    X, Y = np.meshgrid(uperp_vA,magnetic_field_angle)
    Z = X*np.cos(Y)
    
    # calculate levels according to Delta array
    if not _levels:
        Delta = [0.5]
    else:
        Delta=np.arange(_Deltaspacing,1,_Deltaspacing)
    levels=np.abs([-i/(l-i) for i in Delta])

    # plot contours
    fig,ax=plt.subplots(figsize=(8,5))
    # im=ax.imshow(Z,extent=extent,**imkwargs)
    ax1=plt.contourf(X,Y,Z,levels=20,cmap='bwr')
    contours = ax.contour(np.abs(Z),levels,extent=extent,colors='k')

    # reassign contour level labels by Delta parameter
    contours.levels = Delta
    ax.clabel(contours, inline=1,fontsize=10) #, manual=locations
    ax.axhline(np.pi/2,color='k',linestyle='--')

    # colorbar
    # fig.colorbar(im)
    plt.colorbar()
    
    fig.supylabel(r'$\theta$',**tnrfont)
    fig.supxlabel(r'$u_\parallel/v_A$',**tnrfont)
    ax.set_title(r'$l={}$'.format(l),**tnrfont)

    # magnetic field angle limits
    # _yticks = [-np.pi,-2*np.pi/3,-np.pi/3,0,np.pi/3,np.pi/2,2*np.pi/3,np.pi]
    # _yticklabels = [r'$-\pi$',r'$-2\pi/3$',r'$-\pi/3$',r'$0$',r'$\pi/3$',r'$\pi/2$',r'$2\pi/3$',r'$\pi$']
    # ax.set_yticks(_yticks)
    # ax.set_yticklabels(_yticklabels)
    # ax.set_ylim(0)   

    fig.savefig('significant_dopplershift_contours_l_{}.png'.format(l),bbox_inches='tight')
    # plt.show()


if __name__=='__main__':
    for l in [1,2,3,4,5,6,7,8,9,10]:
        plot_significant_dop(l=l,rowlim=(np.pi/2-0.1,np.pi/2+0.1))#np.pi/3,2*np.pi/3))
