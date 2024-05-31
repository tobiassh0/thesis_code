
from func_load import *
from matplotlib.collections import LineCollection
import matplotlib.colors as colors
import random

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

def tracer_plotPhaseSpace(loc,species=['Deuterons'],colors=['k'],num_particles=10,plot=False,phasespace=False):
    simloc=getSimulation(loc)
    ind_lst=list_sdf(simloc)
    xid = np.zeros((len(species),num_particles,len(ind_lst)))
    vxid = np.zeros((len(species),num_particles,len(ind_lst)))
    times = np.zeros(len(ind_lst))
    IDarr = []
    for k in range(len(species)):
        for j in range(num_particles):
            # load files
            for i in range(0,len(ind_lst)):
                d=sdfread(i)
                times[i]+=d.__dict__['Header']['time']
                # for k in getKeys(d): print(k)
                id = getQuantity1d(d,'Particles_ID_subset_tracer_'+species[k])
                if i == 0:
                    ID = random.choice(id)
                    IDarr.append(ID)
                    print(ID)
                index = np.where(id==ID)[0][0]
                vx = d.__dict__['Particles_Px_subset_tracer_'+species[k]]
                x = vx.grid.data[0]
                xid[k,j,i] = x[index]
                vxid[k,j,i]= vx.data[index]
                if plot:
                    if phasespace:
                        # phase-space scatter
                        plt.scatter(xid[k,j,i],vxid[k,j,i],color=colors[k],alpha=(i/len(ind_lst))**4) # deuterons (blue)
            if plot:
                if phasespace:
                    # phase-space lines
                    plt.plot(xid[k,j,:],vxid[k,j,:],'-k',alpha=0.5) # color
                else:
                    # velocity through time
                    plt.plot(ind_lst,vxid[k,:,:],color=colors[k])
    return xid, vxid, times/(num_particles*len(species)), IDarr

if __name__=='__main__': 
    # # test two stream case
    # tracer_twoStream_test(loc='/home/space/phrmsf/Documents/EPOCH/epoch_ID/epoch1d/tracer')

    # D_He3 5% tracer test case
    home = '/home/space/phrmsf/Documents/EPOCH/epoch_ID/epoch1d/D_He3_tracer'
    species = ['Deuterons','He3']
    colors = ['b','r']
    xid,vxid,times,_=tracer_plotPhaseSpace(loc=home,species=species,colors=colors,plot=True,phasespace=True,num_particles=5)
    plt.show()
    print(xid.shape,vxid.shape)
    wcp=(const.qe*3.7/getMass('Protons'))
    tcp = 2*const.PI/wcp
    # # plot first species, first tracked particle, all time
    # plt.plot(times/tcp,xid[0,0,:])
    # plt.show()
    # # plot first species, all tracked particles, all time
    # for i in range(xid.shape[1]):
    #     plt.plot(times/tcp,xid[0,i,:]-np.mean(xid[0,i,:10]))
    # plt.show()
    # # summated FFT spectra on all particles
    # timestep = (times[-1]-times[0])/len(times)
    # n=vxid[0,0,:].size
    # fnyq=0.5*2*const.PI/timestep
    # freq=np.linspace(-fnyq,fnyq,n) 
    # FFT=0
    # for i in range(vxid.shape[1]):
    #     print(i)
    #     pre=np.fft.fft(vxid[0,i,:])
    #     FFT+=np.fft.fftshift(pre)
    # plt.plot(freq/wcp,np.log10(np.abs(FFT))) ; plt.show()