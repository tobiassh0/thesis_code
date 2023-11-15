

## Calculates the phase correlation between two signals, sig & sig0 (only need fft)
# returns the value of the shift corresponding to the phase difference between both arrays
def phaseCorrelation(sig,fft_sig0,dt,tnorm,tmax=35):
	fft_sig = np.fft.fft(sig)
	fft_sig_conj = np.conj(fft_sig)
	R = (fft_sig0 * fft_sig_conj) / abs(fft_sig0 * fft_sig_conj)
	r = np.fft.ifft(R)
	shift = dt*np.argmax(r)
	if shift > dt*len(r)/2:
		shift = shift-tmax*tnorm
	return shift

if __name__=='__main__':
	import numpy as np
	x=np.linspace(0,2*const.PI,100)
	sig0=np.cos(x) ; sig=np.sin(x)
	phaseCorrelation(sig,np.fft(sig0),x[-1]/len(x),1,tmax=2*const.PI)

#--------------------------------------------------------------------#
#
### phase correlation example
#t = np.linspace(0,10,1000)
#f = 2*const.PI
#
#noise_frac = 0.75
#
### signals
#a1 = 1. ; b1 = 1
#a2 = 6. ; b2 = b1
#noise_floor = noise_frac*b1
#pa = 0 ; ps = 0
#noise1 = 0 # np.random.normal(0,noise_floor,len(t)) + 5
#noise2 = 0 # np.random.normal(0,noise_floor,len(t)) + 5
#sig1 = np.exp(-b1*(t-a1)**2)+noise1
#sig2 = np.exp(-b2*(t-a2)**2)+noise2
##sig1 = noise1 ; sig2 = noise2
#
### phase correlation
#shift = phaseCorrelation(sig1,fft_sig2,dt=t[-1]/len(t),tnorm=1,tmax=t[-1])
#plt.plot(t,shift) ; plt.show()

# ------------------------------------ #
### further example
#fft_sig1 = np.fft.fft(sig1)
#fft_sig2 = np.fft.fft(sig2)
#fft_sig1_conj = np.conj(fft_sig1)
#freqs = np.fft.fftfreq(t.shape[-1])
#R = (fft_sig2 * fft_sig1_conj) / abs(fft_sig1 * fft_sig1_conj)
#r = np.fft.ifft(R)
#fig,axs=plt.subplots(nrows=2,sharex=True)
#fig.subplots_adjust(hspace=0.15)
#axs[1].ticklabel_format(useOffset=False)
#axs[0].plot(t,sig1,'r',t,sig2,'b')
#axs[1].plot(t[1:],r[1:],'k')
#
### vertical lines
#ymin = -0.0012
#rmax = np.argmax(r)
#rshift = np.abs(r[rmax])
#tshift = t[rmax]
#print(tshift,rshift)
#axs[1].scatter(tshift,rshift,facecolor='k')
#axs[1].plot([tshift,tshift],[ymin,rshift],linestyle='--',color='darkgrey')
#
### spread second-signal
#b2 = 5
#sig2 = np.exp(-b2*(t-a2)**2)
#fft_sig1 = np.fft.fft(sig1)
#fft_sig2 = np.fft.fft(sig2)
#fft_sig1_conj = np.conj(fft_sig1)
#freqs = np.fft.fftfreq(t.shape[-1])
#R = (fft_sig2 * fft_sig1_conj) / abs(fft_sig1 * fft_sig1_conj)
#r = np.fft.ifft(R)
#
#axs[0,1].plot(t,sig1,'r',t,sig2,'b')
#axs[1,1].plot(t[1:],r[1:])
#
## formatting and labelling
#axs[1].set_ylim(ymin,0.003)
##axs[1,1].set_ylim(-13,13)
#axs[1].set_xlabel(r'$x$',**tnrfont)
#axs[0].set_ylabel(r'Signals',**tnrfont)
#axs[1].set_ylabel(r'Phase-correlation',**tnrfont)
##fig.savefig('thesis/phase-correlation.png')
#plt.show()

