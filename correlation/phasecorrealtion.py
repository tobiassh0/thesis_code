from list_new import *

def phaseCorrelation(sig,fft_sig0,dw,wnorm,wmax=35):
	fft_sig = np.fft.fft(sig)
	fft_sig_conj = np.conj(fft_sig)
	R = (fft_sig0 * fft_sig_conj) / abs(fft_sig0 * fft_sig_conj)
	r = np.fft.ifft(R)
	shift = dw*np.argmax(r)
	if shift > dw*len(r)/2:
		shift = shift-wmax*wnorm
	return shift

#t=[]
#tf=[]
#tf.append(0)
#xiT = np.array([0,0.01,0.05,0.11,0.18,0.3,0.5])
#labels = ['D_99_T_01','D_95_T_05','D_89_T_11','D_82_T_18','D_70_T_30','0_50']
#sim0 = getSimulation('/storage/space2/phrmsf/traceT_0_00')
#pow0 = 10**read_pkl('log10_power')
#fft_pow0 = np.fft.fft(pow0)
#fft_pow0_conj = np.conj(fft_pow0)
#R0 = (fft_pow0 * fft_pow0_conj) / abs(fft_pow0 * fft_pow0_conj)
#r0 = np.fft.ifft(R0)
#shift0 = np.argmax(r0)
#t.append(shift0)
#height = 9 ; width = 7
#fig, ax = plt.subplots(nrows=2,figsize=(width,height))#,sharex=True)
#fig.subplots_adjust(hspace=0.3)
#wcD = const.qe*2.1/getMass('Deuterons')
#for T3 in labels:
#	simloc = getSimulation('/storage/space2/phrmsf/traceT_'+T3)
#	powT = 10**read_pkl('log10_power')
##	wcD = getCyclotronFreq(sdfread(0),'Deuterons')
#	w = read_pkl('omegas_power')
#	dw = ((w[-1]-w[0])/len(w))

#	fft_powT = np.fft.fft(powT)
#	fft_powT_conj = np.conj(fft_powT)

#	R = (fft_pow0 * fft_powT_conj) / abs(fft_pow0 * fft_powT_conj)
#	r = np.fft.ifft(R)
#	freqs = dw*np.arange(0,len(r),1)
#	arg = np.argmax(r)
#	shift = freqs[arg]

#	ax[0].plot(-freqs/wcD,r)
##	plt.scatter(freqs[shift],r[shift])
#	if shift > freqs[-1]/2:
#		shift = shift-freqs[-1]
#	t.append(-shift)
#	tshift = phaseCorrelation(powT,fft_pow0,dw,wcD)
#	tf.append(-tshift)
##	print('time shift = %d' % (shift))

#t = np.array(t) ; tf = np.array(tf)
#fit = np.polyfit(xiT,t,deg=1)
#ax[1].plot(xiT,(fit[0]*xiT+fit[1])/wcD,linestyle='--',color='red')
#ax[1].scatter(xiT,t/wcD,color='k')
##ax[1].scatter(xiT,tf/wcD,marker='x')
#print(t/wcD,tf/wcD)

#ax[0].set_xlabel(r'$\omega_{off}/\Omega_D$',fontsize=20)
#ax[1].set_xlabel(r'$\xi_T$',fontsize=20)
#ax[0].set_ylabel(r'Phase Correlation',fontsize=20)
#ax[1].set_ylabel(r'$\omega_{off}/\Omega_D$',fontsize=20)
#ax[1].set_xlim(0,0.5)# ; ax[1].set_ylim(
#print(fit/wcD)
#fig.savefig('/storage/space2/phrmsf/paper/power-comparison/power_phasecorr.png',bbox_inches='tight')
#fig.savefig('/storage/space2/phrmsf/paper/power-comparison/power_phasecorr.eps',bbox_inches='tight')
#plt.show()

#--------------------------------------------------------------------#

#--------------------------------------------------------------------#

## phase correlation example
t = np.linspace(0,10,1000)
f = 2*const.PI

noise_frac = 0.75

## signals
a1 = 1. ; b1 = 1
a2 = 6. ; b2 = b1
noise_floor = noise_frac*b1
pa = 0 ; ps = 0
noise1 = 0 # np.random.normal(0,noise_floor,len(t)) + 5
noise2 = 0 # np.random.normal(0,noise_floor,len(t)) + 5
sig1 = np.exp(-b1*(t-a1)**2)+noise1
sig2 = np.exp(-b2*(t-a2)**2)+noise2
#sig1 = noise1 ; sig2 = noise2

## phase correlation
fft_sig1 = np.fft.fft(sig1)
fft_sig2 = np.fft.fft(sig2)
fft_sig1_conj = np.conj(fft_sig1)
freqs = np.fft.fftfreq(t.shape[-1])
R = (fft_sig2 * fft_sig1_conj) / abs(fft_sig1 * fft_sig1_conj)
r = np.fft.ifft(R)
fig,axs=plt.subplots(nrows=2,sharex=True)
fig.subplots_adjust(hspace=0.15)
axs[1].ticklabel_format(useOffset=False)
axs[0].plot(t,sig1,'r',t,sig2,'b')
axs[1].plot(t[1:],r[1:],'k')

## vertical lines
ymin = -0.0012
rmax = np.argmax(r)
rshift = np.abs(r[rmax])
tshift = t[rmax]
print(tshift,rshift)
axs[1].scatter(tshift,rshift,facecolor='k')
axs[1].plot([tshift,tshift],[ymin,rshift],linestyle='--',color='darkgrey')

## spread second-signal
b2 = 5
sig2 = np.exp(-b2*(t-a2)**2)
fft_sig1 = np.fft.fft(sig1)
fft_sig2 = np.fft.fft(sig2)
fft_sig1_conj = np.conj(fft_sig1)
freqs = np.fft.fftfreq(t.shape[-1])
R = (fft_sig2 * fft_sig1_conj) / abs(fft_sig1 * fft_sig1_conj)
r = np.fft.ifft(R)

axs[0,1].plot(t,sig1,'r',t,sig2,'b')
axs[1,1].plot(t[1:],r[1:])

# formatting and labelling
axs[1].set_ylim(ymin,0.003)
#axs[1,1].set_ylim(-13,13)
axs[1].set_xlabel(r'$x$',**tnrfont)
axs[0].set_ylabel(r'Signals',**tnrfont)
axs[1].set_ylabel(r'Phase-correlation',**tnrfont)
#fig.savefig('thesis/phase-correlation.png')
plt.show()

#--------------------------------------------------------------------#

#--------------------------------------------------------------------#
simloc = getSimulation('/storage/space2/phrmsf/D_He3_0_10_min_p_0_9') # p_B11_01
Bz = load_batch_fieldmatrix([],'Magnetic_Field_Bz',para=False)
Ex = load_batch_fieldmatrix([],'Electric_Field_Ex',para=False)

L = getGridlen(sdfread(0))
dx = getdxyz(sdfread(0))
t = read_pkl('times')
T = t[-1]
species = getIonSpecies(sdfread(0))
tcmin = 2*const.PI/getCyclotronFreq(sdfread(0),species[-1])
wcmaj = getCyclotronFreq(sdfread(0),species[0])

vA = getAlfvenVel(sdfread(0))
kmax = 0.5*2*const.PI/dx
FT1d = read_pkl('FT_1d_Magnetic_Field_Bz')
print(FT1d.shape)
karg = np.argmax(FT1d[1:,:],axis=1)
karr = np.linspace(0,kmax,FT1d.shape[1])
k_star = karr[karg]
thresh = k_star < 200*wcmaj/vA
k_star = stats.mode(k_star[thresh])[0][0]
dphase = np.linspace(-2*const.PI,2*const.PI,500)
Dx = dphase/k_star # shift between +-2pi, in terms of physical spacing

Rtdx = np.zeros((Bz.shape[0],len(Dx))) # time-space
for i in range(len(Dx)):
	print('roll :: ',i)
	Bzroll = np.roll(Bz,int(Dx[i]/dx),axis=1)
	Rtdx[:,i] = 0.5*(np.sqrt(const.e0/const.mu0))*np.sum(Bzroll*Ex*dx,axis=1)
dumpfiles(Rtdx,'Rtdx_v2')

Rtdx = read_pkl('Rtdx_v2')
print(Rtdx.shape)
plt.axvline(-const.PI,color='darkgrey',linestyle='--')
plt.axvline(const.PI,color='darkgrey',linestyle='--')
plt.imshow((Rtdx),**kwargs,cmap='jet',extent=[dphase[0],dphase[-1],0,T/tcmin]) ; plt.show()











