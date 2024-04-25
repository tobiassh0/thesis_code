
import os,sys
import numpy as np
import matplotlib.pyplot as plt
from func_load import *

# get the ratio between Ex and Delta Bz field energy densities through time and plot each time trace in a panel 
def Ex_Bz(sims,sim0=[],minspec='Protons',mean_to=10,home=None):
    if not home:
        home = os.getcwd()
    if sim0 == []:
        sim0 = sims[0]
    _=getSimulation(sim0)
    times = read_pkl('times')
    tcmin = 2*const.PI/getCyclotronFreq(sdfread(0),minspec)

    # loop through sims and plot ratio
    fratio = np.zeros((len(sims),len(times)))
    rxiHe3  = np.zeros((len(sims),len(times)))
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
        rxiHe3[i,:] = [xiHe3[i] for r in range(len(times))]
        ax[i].plot(times/tcmin,fratio[i,:],color='k')
        ax[i].axhline(1,color='r',linestyle='--')
        ax[i].annotate(xiHe3[i],xy=(0.95,0.95),xycoords='axes fraction',va='top',ha='right')
    ax[0].set_ylim(0,2.5)
    ax[0].locator_params(axis='y',nbins=5)
    ax[0].locator_params(axis='x',nbins=6)
    ax[0].set_xlim(0,times[-1]/tcmin)
    fig.supxlabel(r'$t/\tau_{cp}$',**tnrfont)
    fig.supylabel(r'$E_x^2/\Delta B_z^2$'+'  '+r'$[\epsilon_0\mu_0]$',**tnrfont)
    fig.savefig(home+'Ex_Bz_time.png',bbox_inches='tight')

    return fratio, rxiHe3

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
    xiHe3 = [float(i.split('_')[1]) for i in sims]

    fratio, rxiHe3 = Ex_Bz(sims,sim0=sim0)
    # plt.scatter(rxiHe3,fratio,alpha=0.1)
    # plt.show()
