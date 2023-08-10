
from func_load import * 

#sim_loc = getSimulation('/home/space/phrmsf/Documents/EPOCH/epoch-4.17.16/epoch1d/old/old_sims/twoStream')
#sim_loc = getSimulation('/storage/space2/phrmsf/old/twoStream')
sim_loc = getSimulation('/home/space/phrmsf/Documents/EPOCH/epoch_older/epoch/epoch1d/tests/twostream')
#for k in getKeys(sdfread(0)):print(k)
#raise SystemExit

ind_lst = list_sdf(sim_loc)
print(ind_lst)
Te = 273
vthE = np.sqrt(2*const.kb*Te/getMass('Electrons'))
dend = sdfread(ind_lst[-1])
tend = dend.__dict__['Header']['time']
tpe = 2*const.PI/getPlasmaFreq(sdfread(0),'Left_Electrons')
print(tend/tpe)
LDe = getDebyeLength(sdfread(0),'Left_Electrons')
L = 5e5/LDe

#fig,axs = plt.subplots(nrows=2,ncols=3,sharex=True,sharey=True)
#axs = axs.ravel()
#fig.subplots_adjust(hspace=0.,wspace=0.)

c=0
tprime = [1,2,3,4,5,6]
tprime = ind_lst
for t in tprime:
	tind = int(len(ind_lst)*(t*tpe)/tend)
	print(tind)
	tfile = sdfread(int(t)) #sdfread(ind_lst[tind])
	vxLeft = getQuantity1d(tfile,'Particles_Vx_Left_Electrons')
	vxRight = getQuantity1d(tfile,'Particles_Vx_Right_Electrons')	
#	vxProtons = getQuantity1d(tfile,'Particles_Vx_Protons')
	
	xLeft = np.linspace(0,L,len(vxLeft))
	xRight = np.linspace(0,L,len(vxRight))
#	xProtons = np.linspace(0,L,len(vxProtons))
#	axs[c].scatter(xLeft,vxLeft/vthE,color='r',s=1)
#	axs[c].scatter(xRight,vxRight/vthE,color='b',s=1)
#	axs[c].scatter(xProtons,vxProtons/vthE,color='g',s=1)
#	c+=1
	plt.scatter(xLeft,vxLeft/vthE,color='r',s=1)
	plt.scatter(xRight,vxRight/vthE,color='b',s=1)
#	plt.scatter(xProtons,vProtons/vthE,color='g',s=1)
	plt.savefig('{}_vx.png'.format(t))
#	plt.show()
	plt.clf()


