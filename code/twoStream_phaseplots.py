
from func_load import * 

#sim_loc = getSimulation('/home/space/phrmsf/Documents/EPOCH/epoch-4.17.16/epoch1d/old/old_sims/twoStream')
#sim_loc = getSimulation('/storage/space2/phrmsf/old/twoStream')
sim_loc = getSimulation('/home/space/phrmsf/Documents/EPOCH/5_devel/epoch1d/2stream')
#for k in getKeys(sdfread(0)):print(k)
#raise SystemExit

ind_lst = list_sdf(sim_loc)
print(ind_lst)
d0 = sdfread(0)
for k in getKeys(d0): print(k)
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
tprime = ind_lst
tfile = []
for t in tprime:
	tind = int(len(ind_lst)*(t*tpe)/tend)
	tfile.append(sdfread(int(t))) #sdfread(ind_lst[tind])

files = tfile
plotPhaseSpace(files,['Left_Electrons','Right_Electrons'])
#makevideo(video_name='vx_x')



