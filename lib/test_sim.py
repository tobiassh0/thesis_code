from func_load import *

testsim = getSimulation('/storage/space2/phrmsf/ECRH_JT60U')
# load files
ind_lst = list_sdf(testsim)
d0 = sdfread(0)
dlast = sdfread(ind_lst[-1])
# keys
keys = getKeys(d0)
#for k in keys: print(k)
# times
files = []
for i in range(len(ind_lst)):
	files.append(sdfread(ind_lst[i]))
times = getTimes(files)
t0 = times[0]
tlast = times[-1]
dt = (tlast-t0)/len(ind_lst)
print(t0,tlast,dt,dt*12000)
## CFL condition
dx = getdxyz(d0)
L = getGridlen(d0)
#CFL = bool(np.sqrt(dx**2)/const.c > dt) ; print(CFL)
# frequencies
wce = getCyclotronFreq(d0,'FElectrons')
wpe = getPlasmaFreq(d0,'FElectrons')
(tce, tpe) = 2*const.PI/np.array([wce, wpe])
print(bool(dt==6*tce/12000))

# vx,y,z ring-beam of minority
species = getAllSpecies(d0)
minspec = species[-1]
vx = getQuantity1d(files[0],'Particles_Px_FElectrons')/const.me
vy = getQuantity1d(files[0],'Particles_Py_FElectrons')/const.me
vz = getQuantity1d(files[0],'Particles_Pz_FElectrons')/const.me
print(len(vx),len(vy),len(vz))
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
frac = 100
Nx = L/dx
C = len(vx)/Nx
print('PPCPS :: ',C)
ax.scatter(vx[::frac]/const.c,vy[::frac]/const.c,vz[::frac]/const.c)
plt.show()


