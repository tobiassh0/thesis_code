import numpy as np
import matplotlib.pyplot as plt 
import my_constants as const
from list_new import *
from scipy import signal
from scipy import stats
from scipy.fft import fftshift
from scipy.optimize import curve_fit
from energy import *
import time

def func_ptheory1(dk,a,b): # fit a 1/deltak relation to power using a,b,c coefficients
	return a/(c+b*dk)

def func_ptheory1(dk,a,b): # fit an exp(-deltak) relation to power using a,b coefficients
	return a*np.exp(-b*dk)

def func_ptheory2(dk,a,b,c): # fit a 1/polynomial deltak relation to power using a,b,c coefficients assuming a set power
	return (1E14)/(a*dk**b + c*dk) 

sim_loc = getSimulation('/storage/space2/phrmsf/traceT_0_00')
times = read_pkl('times')
Ex = getQuantity1d(sdfread(0),'Electric_Field_Ex')
nt, nx = len(times), len(Ex)
_, _, _, Smag = getPoynting(times,nt,nx,min_species='Alphas',plot=True)

############################
#sim_loc = getSimulation('/storage/space2/phrmsf/traceT_0_11')
#FT_2d = read_pkl('FT_2d_Magnetic_Field_Bz')
#times= read_pkl('times')
#(nw, nk) = FT_2d.shape
#print(nw,nk)
#nx, nt = 2*nk, len(times)
#d0 = sdfread(0)
#T = times[-1]
#L = getGridlen(d0)
#### get consts
#species = getIonSpecies(d0)
#wc_maj = getCyclotronFreq(d0,species[0])
#wce = getPlasmaFreq(d0,'Electrons')
#wpi = getPlasmaFreq(d0, species[0])
#wpe = getPlasmaFreq(d0, 'Electrons')
#lambdaD = getDebyeLength(d0,'Electrons')
#va = getAlfvenVel(d0)
#wnorm = wc_maj
#knorm = wnorm/va
#dx = getdxyz(sdfread(0))
#print('dx : {}'.format(dx))
#### normalise lims
#klim, wlim = 0.5*nx*2*const.PI/L, 0.5*nt*2*const.PI/T 
#klim_prime, wlim_prime = klim/knorm, wlim/wnorm
#in_klimprime = 100
#in_wlimprime = 50
#### cut FT_2d to area of interest
#w_lim, k_lim = FT_2d.shape[0]*(in_wlimprime/wlim_prime), FT_2d.shape[1]*(in_klimprime/klim_prime)
#FT_2d = FT_2d[:int(w_lim),:int(k_lim)]
#fig, ax = plot2dTransform(FT_2d, [va,True], in_klimprime, in_wlimprime, Omega_label=r'$\Omega_D$' ,cbar=False, clim=(None,None), cmap='magma')
#### plot horizontal harmonics
##for i in range(0,10):
##	plt.axhline(i,alpha=0.5,color='k',linestyle='--')
#### Cold plasma dispersion
#omegas = wnorm*np.linspace(0,in_wlimprime,10000) # do this so that the range of k near w=0 does not diverge, then shift it down
##ax.plot(omegas/wnorm,omegas/wnorm)# w = vA * k
#k1,k2,k3=coldplasmadispersion(d0,'Deuterons','Electrons',z1=1,z2=1,omegas=omegas) # two solutions to the cold plasma dispersion
#knorm = (va/wc_maj)
#k1,k2,k3 = k1*knorm, k2*knorm, k3*knorm
#thresh = k2 > 0 # threshold the array so it only plots the FAW and not the horizontal line to the 0 parts of the dispersion
#ax.plot(k2[thresh],omegas[thresh]/wnorm,color='k',linestyle='-',alpha=0.75)
#print('k1 : ',k1,'\n','k2 : ',k2)
#ax.set_ylim(0,in_wlimprime)
#ax.set_xlim(0,in_klimprime)

### Cold wave modes
#ax = ColdWaveModes(ax,[wpi, wpe, wc_maj, wce],va,wnorm,LHact=True) ## check list_new to see which modes we can plot (bool)

### w1 & w2 : Buchsbaum 1960, 2 ion resonance
##w1 = np.sqrt(wce**2 + wpe**2)
##ax.axhline(w1/wnorm,color='k',linestyle=':')
##print(w1/wnorm)
#w2 = np.sqrt(wce*wc_maj*((wce*wc_maj + wpe**2)/((wce**2)+(wpe**2))))
#ax.axhline(w2/wnorm,color='k',linestyle='-.')
#print(w2/wnorm)
#plotting(fig,ax,'FT_2d_Magnetic_Field_Bz')

############################
#sim_loc = getSimulation('/storage/space2/phrmsf/cold_OLD')
#quant='Magnetic_Field_Bz'
#fm = read_pkl('fieldmatrix_'+quant)
#times=read_pkl('times')
#(nt,nx) = fm.shape
#karea = 40
#warea = 50
#nfft = nt//10
#noverlap = nfft//2
#T=times[-1]
#dt=T/nt
#L=getGridlen(sdfread(0))
#wc_maj = getCyclotronFreq(sdfread(0),'Deuterons')
#va = getAlfvenVel(sdfread(0))
##wc_maj=100599816.674923
##va=10218477.723612128
##wlim:86817299821.74399
##klim:140894.9327307366
##wlim_prime:862.9966007024082
##klim_prime:14311.474703092077
##wci:100599816.674923
##wpi:2926674062.4692698
##wce:369352226921.9798
##wpe:178398640789.3469

##bispec = read_pkl('Bicohmat_ka_wa_80_40_nfft_2000_noverlap_1000')
##print('Bicoh shape ::',bispec.shape)
##print(bispec)
###min(bispec)
###max(bispec)
##bispec = np.log10(bispec)
##extent=[0,karea,0,karea]
##fig,ax,_=plot_bicoh(bispec,extent=extent,smooth=True,cbar=True,cmap='jet')
##ax.set_xlabel(r'$k_1v_A/\Omega_{D}$',fontsize=18)
##ax.set_ylabel(r'$k_2v_A/\Omega_{D}$',fontsize=18)
##plotting(fig,ax,'bicoh')

#getBicoh(karea,warea,fm,dt,T,L,wc_maj,va,nfft=nfft,noverlap=noverlap,window=True,bispectrum=False,cmap='jet')
#plt.show()

##################################
#from scipy.fftpack import next_fast_len
#from polycoherence import _plot_signal, polycoherence, plot_polycoherence

### signal
#N = 10001
#t = np.linspace(0, 100, N)
#fs = 1 / (t[1] - t[0])
#s1 = np.cos(2 * pi * 5 * t + 0.2)
#s2 = 3 * np.cos(2 * pi * 7 * t + 0.5)
#np.random.seed(0)
#noise = 5 * np.random.normal(0, 1, N)
#signal = s1 + s2 + 0.5 * s1 * s2 + noise
#_plot_signal(t, signal)

#sim_loc = getSimulation('/storage/space2/phrmsf/traceT_highres')
#quant = 'Magnetic_Field_Bz'
#fm = read_pkl('fieldmatrix_'+quant)
#times = read_pkl('times')
#(nt,nx) = fm.shape
#N = nt
#wcyc = getCyclotronFreq(sdfread(0),'Deuterons')
#tc_D = 2*const.PI/wcyc
#tn = 3.5 # normalised
#tind = N*(tn/(times[-1]/tc_D))
#fm = fm[int(tind),:] # all space, one time
#L = getGridlen(sdfread(0))
#x = np.linspace(0,L,nx)
#_plot_signal(x,fm)

#va = getAlfvenVel(sdfread(0))
#fs = 1/getdxyz(sdfread(0))
#knorm = wcyc/va

#flim1 = (0,80*knorm) ; flim2 = flim1
#kw = dict(nperseg=N//5, noverlap=N//10, nfft=next_fast_len(N//2))
#freq1, freq2, bicoh = polycoherence(fm, fs, **kw)
#plot_polycoherence(freq1, freq2, bicoh)

#### calc and plot bicoh
#kw = dict(nperseg=N // 10, noverlap=N // 20, nfft=next_fast_len(N // 2))
#freq1, freq2, bicoh = polycoherence(signal, fs, **kw)
#plot_polycoherence(freq1, freq2, bicoh)

#### calc and plot bispec
#freq1, fre2, bispec = polycoherence(signal, fs, norm=None, **kw)
#plot_polycoherence(freq1, fre2, bispec)

#### calc and plot part of bicoh
#flim1 = (0,10) ; flim2 = flim1
#freq1, freq2, bicoh = polycoherence(signal, fs, flim1=flim1, flim2=flim2, **kw)
#plot_polycoherence(freq1, freq2, bicoh, flim1, flim2)

##################################

### calc & plot bicoherence
#sim_loc = getSimulation('/storage/space2/phrmsf/traceT_noT')
#karea = 80
#warea = 40
##nfft = Nt//5
##noverlap = nfft//2
##bicname = 'Bispecmat_ka_wa_80_40_nfft_2400_noverlap_1200'
#bicname = 'Bicohmat_ka_wa_80_40_nfft_2400_noverlap_1200'
#bic = read_pkl(bicname)
#fig,ax,_=plot_bicoh(bic,extent=[0,karea,0,karea],bispectrum=False,smooth=True,cbar=True,clim=(0,1),cmap='jet')
#ax.set_xlabel(r'$k_1v_A/\Omega_D$',fontsize=18)
#ax.set_ylabel(r'$k_2v_A/\Omega_D$',fontsize=18)

#plotting(fig,ax,'bicoherence_karea_{}_warea_{}'.format(karea,warea))


##################################










