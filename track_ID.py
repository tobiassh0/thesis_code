
from func_load import *
from matplotlib.collections import LineCollection
import matplotlib.colors as colors
import energy.LarmorRadii as lr
import collections
import itertools


# makes dict of allids (shuffled) and a linear index (0-->len(px))
# finds linear index of tracked ids in allids using dict key
def track_IDS(linear_index,allids,trackids=[None]):
    result = dict(zip(allids,linear_index))
    return np.array([result[idx] for idx in trackids])

def tracer_twoStream_test(loc):
    simloc=getSimulation(loc)
    ind_lst = list_sdf(simloc)
    # ID of "Left" particle you want to track
    ID = 2620 # 2631
    # load files
    xid = np.zeros(len(ind_lst))
    vid = np.zeros(len(ind_lst))
    for i in ind_lst:
        d=sdfread(i)
        idLeft = getQuantity1d(d,'Particles_ID_subset_tracer_Left')
        # idRight = getQuantity1d(d,'Particles_ID_subset_tracer_Right')
        index = np.where(idLeft==ID)[0][0]
        # get all velocities
        vxLeft = d.__dict__['Particles_Px_subset_tracer_Left']
        xLeft = vxLeft.grid.data[0]
        vxRight = d.__dict__['Particles_Px_subset_tracer_Right']
        xRight = vxRight.grid.data[0]
        xid[i] = xLeft[index]
        vid[i] = vxLeft.data[index]
        # # plot all particles
        # plt.scatter(xLeft,vxLeft.data,color='b',alpha=0.1)
        # plt.scatter(xRight,vxRight.data,color='r',alpha=0.1)
        # # plot traced particle
        # plt.scatter(xLeft[index],vxLeft.data[index],color='k')
        # # title and save
        # plt.title(i)
        # plt.savefig('trace_particles_{}.png'.format(i))
        # plt.clf()
        plt.scatter(xid[i],vid[i],c='k',s=50,edgecolor='none',alpha=i/len(ind_lst))
    plt.plot(xid,vid,'--b',alpha=0.2)
    plt.xlabel('x') ; plt.ylabel('vx')
    plt.title('Left ID : {}'.format(ID))
    plt.savefig('trace_particle_ID_{}.png'.format(ID))
    return None

def tracer_getXPX(loc,species=['Deuterons'],minspec=['Protons'],colors=['k'],num_particles=10,plot=False,phasespace=False):
    # load sim & set constants
    sim=getSimulation(loc)
    d0=sdfread(0)
    tcmin=2*const.PI/getCyclotronFreq(d0,'Protons')
    # find files where ID is output
    restart_files=lr.para_check_restart(sim) # list_sdf(sim)
    # setup empty arrays for positions, momentum & times
    xid = np.zeros((len(species),num_particles,len(restart_files)))
    pxid = np.zeros((len(species),num_particles,len(restart_files)))
    times = np.zeros(len(restart_files))
    # setup which particles to initially track by randomly choosing
    IDarr = np.zeros((len(species),num_particles)) # array of IDs shape [No. species, No. particles to follow]
    for k in range(len(species)):
        tid = getQuantity1d(d0,'Particles_ID_'+species[k])
        IDarr[k,:] = np.random.choice(tid,num_particles,replace=False)
    # loop through each time, find particle and append array with position and x-momentum (TODO; extend to py & pz)
    for i in range(len(restart_files)):
        print("{:.4f}%".format(100*i/len(restart_files)))
        di=sdfread(restart_files[i])
        times[i]=di.__dict__['Header']['time']
        for k in range(len(species)):
            px=di.__dict__['Particles_Px_'+species[k]]
            x=px.grid.data[0]
            ids=getQuantity1d(di,'Particles_ID_'+species[k])
            linear_index=np.arange(0,len(ids),1)
            indices=track_IDS(linear_index=linear_index,allids=ids,trackids=IDarr[k,:]) # shape: [No. particles]
            for j in range(num_particles): # len(indices)
                index=indices[j]
                pxid[k,j,i]=px.data[index]
                xid[k,j,i]=x[index]
    # plotting scripts
    if plot:
        if phasespace:
            for k in range(len(species)):
                for j in range(num_particles):
                    for i in range(len(times)):
                        plt.scatter(xid[k,j,i],pxid[k,j,i],color=colors[k],alpha=i/len(times)**4)
            plt.xlabel(r'$x$'+'  '+r'$[m]$',**tnrfont)
            plt.ylabel(r'$p_x$'+'  '+r'$[kgms^{-2}]$',**tnrfont)
            plt.legend(species,loc='best')
            plt.savefig('track_phase_space.png')
        else: # velocity space 
            for k in range(len(species)):
                for j in range(num_particles):
                    plt.plot(times/tcmin,pxid[k,j,:],color=colors[k])
            plt.xlabel(r'$t/\tau_{c\sigma}$',**tnrfont)
            plt.ylabel(r'$p_x$'+'  '+r'$[kgms^{-2}]$',**tnrfont)
            plt.legend(species,loc='best')
            plt.savefig('track_momentum_time.png')

    return times, xid, pxid, IDarr


def tracer_plot3dPhaseSpace(loc,species=['Deuterons','He3'],colors=['b','r'],num_particles=10):
    sim=getSimulation(loc)
    restart_files=lr.para_check_restart(sim)
    restart_times=np.zeros(len(restart_files))
    d0=sdfread(0)
    tcp=2*const.PI/getCyclotronFreq(d0,'Protons')
    # two (species) 2d arrays [No. particles (IDs), times, [x,px]]
    ph1=np.zeros((num_particles,len(restart_files),2))
    ph2=np.zeros((num_particles,len(restart_files),2))
    # choose random IDs for each species to follow
    id1 = getQuantity1d(d0,'Particles_ID_subset_tracer_'+species[0])
    id2 = getQuantity1d(d0,'Particles_ID_subset_tracer_'+species[1])
    # normdist to weight id's so that choose particles close to centre of sample
    normdist1 = (1/(3000*(2*const.PI)**0.5))*np.exp(-0.5*((np.linspace(0,len(id1),len(id1))-len(id1)/2)/3000)**2)
    normdist1/=np.sum(normdist1*1) # spacing of 1
    normdist2 = (1/(3000*(2*const.PI)**0.5))*np.exp(-0.5*((np.linspace(0,len(id2),len(id2))-len(id2)/2)/3000)**2)
    normdist2/=np.sum(normdist2*1) # spacing of 1
    # select IDs randomly in list without repeats (replace=F)
    ID1 = np.random.choice(id1,num_particles,replace=False,p=normdist1)
    ID2 = np.random.choice(id2,num_particles,replace=False,p=normdist2)
    for i in range(len(restart_files)):
        d=sdfread(i)
        restart_times[i]=d.__dict__['Header']['time']
        id1 = getQuantity1d(d,'Particles_ID_subset_tracer_'+species[0])
        id2 = getQuantity1d(d,'Particles_ID_subset_tracer_'+species[1])
        for j in range(num_particles):
            index1 = np.where(id1==ID1[j])[0][0]
            index2 = np.where(id2==ID2[j])[0][0]
            px1 = d.__dict__['Particles_Px_subset_tracer_'+species[0]]
            px2 = d.__dict__['Particles_Px_subset_tracer_'+species[1]]
            x1 = px1.grid.data[0]
            x2 = px2.grid.data[0]
            ph1[j,i,:]=[x1[index1],px1.data[index1]]
            ph2[j,i,:]=[x2[index2],px2.data[index2]]
    fig=plt.figure()
    ax=fig.add_subplot(projection='3d')
    for j in range(num_particles):
        x1 = ph1[j,:,0]
        px1 = ph1[j,:,1]
        ax.plot(x1,restart_times/tcp,px1,color=colors[0])
        x2 = ph2[j,:,0]
        px2 = ph2[j,:,1]
        ax.plot(x2,restart_times/tcp,px2,color=colors[1])
    ax.set_xlabel(r'$x$'+'  '+r'$[m]$',**tnrfont)
    ax.set_ylabel(r'$t/\tau_{cp}$',**tnrfont)
    ax.set_zlabel(r'$p_x$'+'  '+r'$[kgm/s^{-2}]$',**tnrfont)
    plt.show()
    fig.savefig('tracer_D_He3_phase_time.png',bbox_inches='tight')
    
    return None


if __name__=='__main__': 
    # # test two stream case
    # tracer_twoStream_test(loc='/home/space/phrmsf/Documents/EPOCH/epoch_ID/epoch1d/tracer')

    # # D_He3 5% tracer test case
    # home = '/home/space/phrmsf/Documents/EPOCH/epoch_ID/epoch1d/D_He3_tracer'
    # species = ['Deuterons']#,'He3']
    # colors = ['b','r']
    # tracer_plot3dPhaseSpace(loc=home,num_particles=20)
    # times,xid,pxid,_=tracer_getXPX(loc=home,species=species,colors=colors,plot=True,phasespace=False,num_particles=1)
    # # plt.show()

    # # D_He3 high res tracers
    home = '/storage/space2/phrmsf/lowres_D_He3/tracer_0_05/'
    species=['Deuterons','He3']
    colors=['b','r']
    num_particles = 1 # particles to track for each species
    times,xid,pxid,_=tracer_getXPX(loc=home,species=species,colors=colors,plot=True,phasespace=False,num_particles=num_particles)
    # Find FT
    dt = (times[-1]-times[0])/len(times)
    wnyq=0.5*2*const.PI/dt
    freq=np.linspace(-wnyq,wnyq,len(times))/getCyclotronFreq(sdfread(0),'Protons')
    # loop through each species and particle id
    fig,ax=plt.subplots(figsize=(10,6),nrows=len(species))
    for k in range(len(species)):
        for j in range(num_particles):
            window=np.hanning(len(pxid[k,j,:]))
            FT=np.fft.fftshift(
                np.fft.fft(window*pxid[k,j,:])
                )
            ax[k].plot(freq,np.abs(FT),color=colors[k])
        ax[k].set_ylabel('FT magnitude',**tnrfont)
        ax[k].set_xlim(0) # set minimum above 0
    fig.supxlabel(r'$\omega/\Omega_p$',**tnrfont)
    fig.savefig('FT_momentum.png',bbox_inches='tight')
    plt.clf()
