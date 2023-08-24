
from func_load import *

sim = 'ECRH/ECRH_JT60U_3'
home = '/storage/space2/phrmsf/'

## setup
fig, ax = plt.subplots(figsize=(8,4)) #(8,4)
sim_loc = getSimulation(home+sim)
d0 = sdfread(0)

FT_2d = read_pkl('FT_2d_Magnetic_Field_Bz')
times = read_pkl('times')
dt = times[-1]/len(times)
dx = getdxyz(d0)
LDe = getDebyeLength(d0,'Electrons')
## freq limits
klim = 2*0.5*const.PI/dx
wlim = 2*0.5*const.PI/dt
vA = getAlfvenVel(d0)
wcyce = getCyclotronFreq(d0,'Electrons')
wcycp= getCyclotronFreq(d0,'Protons')
wnorm = wcycp
knorm = wcycp/vA
klim_prime = klim/knorm
wlim_prime = wlim/wnorm
print(klim_prime,wlim_prime)
(nw,nk) = FT_2d.shape
print(nw,nk)
#	kmax = 40; wmax = 100
kmax = klim_prime; wmax = wlim_prime
dk = 2*klim/(2*nk) ; dw = (2*wlim/len(times))

#	_,_ = power(klim_prime=klim_prime,wlim_prime=wlim_prime,wmax=35,kmax=kmax,quantity='Magnetic_Field_DeltaBz',plot=False,read=False)

# print dw, dk, klim and wlim in various units
print('klim (vA/omega_p, va/omega_e) :: ',klim*vA/wcycp,klim*vA/wcyce) 
print('wlim (omega_p, omega_e) :: ',wlim/wcycp,wlim/wcyce)
print('dk (vA/omega_p, va/omega_e) :: ',dk*vA/wcycp,dk*vA/wcyce) 
print('dw (omega_p, omega_e) :: ',dw/wcycp,dw/wcyce) 	

# chopping and plotting
FT_2d = FT_2d[:int(nw*wmax/wlim_prime),:int(nk*kmax/klim_prime)]
ax.imshow(np.log10(FT_2d),interpolation='nearest',cmap='magma',origin='lower',aspect='auto',extent=[0,kmax,0,wmax],vmin=-4,vmax=6)
ax.set_xlabel(r'$kv_A/\Omega_e$',fontsize=18) # v_A/\Omega_D
ax.set_ylabel(r'$\omega/\Omega_e$',fontsize=18)
plt.show()
