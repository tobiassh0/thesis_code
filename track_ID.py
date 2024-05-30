
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

if __name__=='__main__': 
    # # test two stream case
    # tracer_twoStream_test(loc='/home/space/phrmsf/Documents/EPOCH/epoch_ID/epoch1d/tracer')
    home = '/home/space/phrmsf/Documents/EPOCH/epoch_ID/epoch1d/D_He3_tracer'
    simloc=getSimulation(home)
    ind_lst = list_sdf(simloc)
    name1 = 'Deuterons' ; name2 = 'He3'
    # ID of "name" particle you want to track
    ID = 35 # default but should choose based on subset
    # load files
    for j in range(10):
        x1id = []
        vx1id = []
        x2id = []
        vx2id = []
        for i in range(0,len(ind_lst)-1):
            d=sdfread(i)
            # for k in getKeys(d): print(k)
            id1 = getQuantity1d(d,'Particles_ID_subset_tracer_'+name1)
            id2 = getQuantity1d(d,'Particles_ID_subset_tracer_'+name2)
            if i == 0:
                ID1 = random.choice(id1)
                ID2 = random.choice(id2)
                print(ID1, ID2)
            index1 = np.where(id1==ID1)[0][0]
            index2 = np.where(id2==ID2)[0][0]
            vx1 = d.__dict__['Particles_Px_subset_tracer_'+name1]
            vx2 = d.__dict__['Particles_Px_subset_tracer_'+name2]
            x1 = vx1.grid.data[0]
            x2 = vx2.grid.data[0]
            x1id.append(x1[index1])
            vx1id.append(vx1.data[index1])
            x2id.append(x2[index2])
            vx2id.append(vx2.data[index2])
            # phase-space scatter
            plt.scatter(x1id[-1],vx1id[-1],color='b',alpha=(i/len(ind_lst))**4) # deuterons (blue)
            plt.scatter(x2id[-1],vx2id[-1],color='r',alpha=(i/len(ind_lst))**4) # helium3 (red)
        # phase-space lines
        plt.plot(x1id,vx1id,'-k',alpha=0.5)
        plt.plot(x2id,vx2id,'-k',alpha=0.5)
        # # velocity through time
        # plt.plot(ind_lst[:-1],vx1id,color='b')
        # plt.plot(ind_lst[:-1],vx2id,color='r')
    plt.show()
