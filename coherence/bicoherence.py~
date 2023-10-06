from func_load import *

simloc=getSimulation('/storage/space2/phrmsf/ECRH/ECRH_JT60U_5_2')
d0 = sdfread(0)

fieldmatrix = load_batch_fieldmatrix([],'Magnetic_Field_Bz')
times = read_pkl('times')
L = getGridlen(d0)
Nt = len(times)
T = times[-1]
dt = T/Nt

wnorm = getCyclotronFreq(d0,'Protons')
knorm = wnorm/getAlfvenVel(d0)
bispec = True
klabel = r'$v_A/\Omega_p$'

nfft = Nt//5
noverlap = nfft//2
karea = 30
warea = 30
fig, ax = getBicoh(karea,warea,fieldmatrix,dt,T,L,wnorm,knorm,nfft=nfft,\
		noverlap=noverlap,window=True,bispectrum=bispec,klabel=klabel)

