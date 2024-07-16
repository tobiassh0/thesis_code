
import os,sys
import numpy as np
import matplotlib.pyplot as plt
from func_load import *

# get the ratio between Ex and Delta Bz field energy densities through time and plot each time trace in a panel 
def Ex_Bz(sims,xarr,sim0=[],minspec='Protons',mean_to=10,home=None):
    """
        In:
            sims        : list of sims to loop over and compare (must have len > 1)
            xarr        : x array which is looped over, same len as sims (i.e. xiHe3)
            minspec     : minority species to normalise time to
            mean_to     : # time files to take mean to
            home        : location to save fig, if None then save in cwd
        Out:
            fratio      : ratio between field energy densities Ex/Bz
            rxarr       : matching coordinates for each ratio through time per xarr[i]
    """
    if not home:
        home = os.getcwd()
    if sim0 == []:
        sim0 = sims[0]
    _=getSimulation(sim0)
    times = read_pkl('times')
    tcmin = 2*const.PI/getCyclotronFreq(sdfread(0),minspec)

    # loop through sims and plot ratio
    fratio = np.zeros((len(sims),len(times)))
    rxarr  = np.zeros((len(sims),len(times)))
    M = 4
    N = int(len(sims)/M)
    fig,ax=plt.subplots(figsize=(10,6),ncols=M,nrows=N,sharex=True,sharey=True,layout='constrained')
    # fig.subplots_adjust(hspace=0,wspace=0)
    ax=ax.ravel()
    for i in range(len(sims)):
        simloc = getSimulation(home+sims[i])
        times= read_pkl('times')
        Ex = read_pkl('Exenergy')/(2/const.e0)
        Bz = read_pkl('Bzenergy')/(2*const.mu0)
        DeltaEx = Ex-np.mean(Ex[:mean_to])
        DeltaBz = Bz-np.mean(Bz[:mean_to])
        # plt.plot(DeltaEx,color='r')
        # plt.plot(DeltaBz,color='b')
        # fratio[i] = np.max(DeltaEx)/np.max(DeltaBz)
        fratio[i,:] = DeltaEx/DeltaBz
        rxarr[i,:] = [xarr[i] for r in range(len(times))]
        ax[i].plot(times/tcmin,fratio[i,:],color='k')
        ax[i].axhline(1,color='r',linestyle='--')
        ax[i].annotate("{:.0f}".format(xarr[i]),xy=(0.95,0.95),xycoords='axes fraction',va='top',ha='right')
    ax[0].set_ylim(0,2.5)
    ax[0].locator_params(axis='y',nbins=5)
    ax[0].locator_params(axis='x',nbins=6)
    ax[0].set_xlim(0,times[-1]/tcmin)
    fig.supxlabel(r'$t\Omega_p/2\pi$',**tnrfont)
    fig.supylabel(r'$(E_x^2/\Delta B_z^2)\epsilon_0\mu_0$',**tnrfont)
    fig.savefig(home+'Ex_Bz_time.png',bbox_inches='tight')

    return fratio, rxarr

def f1_f2(home,sims,xarr,f1='Electric_Field_Ex',f2='Magnetic_Field_Bz',sim0=[],minspec='Protons',mean_to=10,ylabel=None,\
            smooth=False,smooth_N=100):
    """
    Finds the ratio between two field energy densities, these are default to Ex and oscillatory Bz. Plots as a multipanel (nrows and
    ncols harcoded for now) for each sim, with approproate label in the top right. Has ability to "smooth" traces through np.convolve.
    Xlabel not given as will return t\Omega_{minspec}/2\pi automatically.
        In:
            home            : str, location to save fig, if None then save in cwd
            sims            : lst, list of sims to loop over and compare (must have len > 1)
            xarr            : np.arr, x array which is looped over, same len as sims (i.e. xiHe3)
            f1              : str, name of first field
            f2              : str, name of second field
            minspec         : str, minority species to normalise time to
            mean_to         : int, number of time files to take mean to
            smooth          : boolean, determines whether want to smooth field energy dens ratios or leave whole
            smooth_N        : int, number used to smooth fratio
        Out:
            fratio          : np.arr, ratio between field energy densities f1/f2 for each sim at all time
            rxarr           : 2d lst, matching coordinates for each ratio through time per xarr[i]
    """
    f1name,f1label,f1mult = Endict(f1)
    f2name,f2label,f2mult = Endict(f2) # name of energy field
    if not home:
        home = os.getcwd()
    if sim0 == []:
        sim0 = sims[0]
    if ylabel == None:
        ylabel = '('+f1label+')'+r'$/$'+'('+f2label+')'

    _=getSimulation(sim0)
    times = read_pkl('times')
    tcmin = 2*const.PI/getCyclotronFreq(sdfread(0),minspec)

    # loop through sims and plot ratio
    fratio = np.zeros((len(sims),len(times)))
    rxarr  = np.zeros((len(sims),len(times)))
    M = 4
    N = int(len(sims)/M)
    fig,ax=plt.subplots(figsize=(10,6),ncols=M,nrows=N,sharex=True,sharey=True,layout='constrained')
    # fig.subplots_adjust(hspace=0,wspace=0)
    ax=ax.ravel()
    for i in range(len(sims)):
        print(sims[i])
        simloc = getSimulation(home+sims[i])
        times= read_pkl('times')
        f1En = read_pkl(f1name)*f1mult
        f2En = read_pkl(f2name)*f2mult
        Deltaf1En = f1En-np.mean(f1En[:mean_to])
        Deltaf2En = f2En-np.mean(f2En[:mean_to])
        # plt.plot(DeltaEx,color='r')
        # plt.plot(DeltaBz,color='b')
        # fratio[i] = np.max(DeltaEx)/np.max(DeltaBz)
        fratio[i,:] = Deltaf1En/Deltaf2En
        rxarr[i,:] = [xarr[i] for r in range(len(times))]
        # smooth traces
        if smooth:
            fieldratio = np.convolve(fratio[i,:],np.ones(smooth_N)/smooth_N,mode='valid')
            ttimes = np.linspace(0,max(times),len(fieldratio))
        else:
            fieldratio = fratio[i,:]
            ttimes = times
        ax[i].plot(ttimes/tcmin,fieldratio,color='k')
        ax[i].axhline(1,color='r',linestyle='--')
        ax[i].annotate(xarr[i],xy=(0.95,0.95),xycoords='axes fraction',va='top',ha='right')
    ax[0].set_ylim(0,3.0)
    ax[0].locator_params(axis='y',nbins=5)
    ax[0].locator_params(axis='x',nbins=6)
    ax[0].set_xlim(0,times[-1]/tcmin)
    fig.supxlabel(r'$t$'+getOmegaLabel(minspec)+r'$/2\pi$',**tnrfont)
    fig.supylabel(ylabel,**tnrfont)
    fig.savefig(home+f1name+'_'+f2name+'_time.png',bbox_inches='tight')
    return fratio, rxarr


if __name__=='__main__':
    from func_load import *
    home='/storage/space2/phrmsf/lowres_D_He3/'
    sims=[]
    for sim in os.listdir(home):
        if 'p_90' in sim:
            sims.append(sim)
    sims = np.sort(sims)
    sim0 = getSimulation(home+sims[0])
    sims = sims[1:]
    xiHe3 = [int(i.split('_')[1]) for i in sims]

    # fratio, rxiHe3 = Ex_Bz(sims,xarr=xiHe3,sim0=sim0,home=home)
    fratio, rxiHe3 = f1_f2(home,sims,xiHe3,'Electric_Field_Ex','Magnetic_Field_Bz',sim0=sim0,\
                            ylabel=r'$(E_x^2/\Delta B_z^2)\epsilon_0\mu_0$',smooth=True)
    # plt.scatter(rxiHe3,fratio,alpha=0.1)
    # plt.show()
