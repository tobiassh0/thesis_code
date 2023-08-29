from func_load import *

simloc = getSimulation('/storage/space2/phrmsf/ECRH/ECRH_JT60U_3')
d0 = sdfread(0)
d1 = sdfread(1)

# check times
dt = getdt(None,d0,d1)
wce = getCyclotronFreq(d0,'Electrons')
wcp = getCyclotronFreq(d0,'Protons')
wpe = getPlasmaFreq(d0,'Electrons')
wpi = getPlasmaFreq(d0,'Protons')
tce = 2*const.PI/wce
tcp = 2*const.PI/wcp
times = read_pkl('times')
print(dt/tce,dt/tcp)

# check freq lim
wUH = getUpperHybrid(wce,wpe)
wlim = 0.5*2*const.PI/dt
dw = 2*wlim/len(times)
print(wlim/wce,wlim/wcp,wlim/wUH)
print(dw/wce,dw/wcp)
