

def shared_area_noise():
	## shared area example
	
	# x-dimension array (arb)
	t = np.linspace(0,10,1000)
	dt = t[-1]/len(t)
	tcc = np.linspace(-t[-1],t[-1],2*len(t)-1)
	f = 2*const.PI
	
	# noise amplitude values to loop through
	vals = np.linspace(0,3.5,50)#[0.,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0,3.25,3.5]
	
	# number of runs to do at each noise amplitude
	N=1e4
	peak_ngauss=np.zeros((len(vals),int(N)))
	peak_gauss=np.zeros((len(vals),int(N)))
	peak_ccorr=np.zeros((len(vals),int(N)))
	peak_ccgauss=np.zeros((len(vals),int(N)))
	
	# signal parameters
	a1 = 1. ; b1 = 1
	a2 = 6. ; b2 = b1
	
	# Monte Carlo method
	for j in range(len(vals)):
		noise_floor = vals[j]*b1
		print(j)
		png = np.zeros(int(N)) ; pg = np.zeros(int(N))
		pcc = np.zeros(int(N)); pccg = np.zeros(int(N))
		for i in range(int(N)):
			# noise and signals	
			noise1 = np.random.normal(0,noise_floor,len(t))# + 5 # creates random noise with each iteration
			noise2 = np.random.normal(0,noise_floor,len(t))# + 5
			sig1 = np.exp(-b1*(t-a1)**2)+noise1
			sig2 = np.exp(-b2*(t-a2)**2)+noise2
	
			# shared area between signals
			area, tpg = shared_area(sig1,sig2,fitgauss=True,dx=dt) # from list_new
			tpng = t[np.argmax(area)]

			# cross correlation
			cc = np.correlate(sig2,sig1,mode='full')
			tpcc = tcc[np.argmax(cc)]
			# # fit Gauss
			# popt, pcov = curve_fit(lambda t,a,b,c: np.exp(a*(t-b)**2)+c,tcc,cc,maxfev=10) # exponential fitting
			# tpccg = popt[1]#TODO;

			# summing total
			png[i]=tpng
			pg[i]=tpg
			pcc[i]=tpcc
			# pccg[i]=tpccg
		
		# append all values in 2d array
		peak_ngauss[j,] = (png)
		peak_gauss[j,] = (pg)
		peak_ccorr[j,] = (pcc)
		# peak_ccgauss[j,] = (pccg)

	# mean and std of each run
	mean_ngauss=np.mean(peak_ngauss,axis=1)
	std_ngauss=np.std(peak_ngauss,axis=1)
	mean_gauss=np.mean(peak_gauss,axis=1)
	std_gauss=np.std(peak_gauss,axis=1)
	mean_cc=np.mean(peak_ccorr,axis=1)
	std_cc=np.std(peak_ccorr,axis=1)
	# mean_ccgauss=np.mean(peak_ccgauss,axis=1)
	# std_ccgauss=np.std(peak_ccgauss,axis=1)
	
	# setup plot
	fig,axs=plt.subplots(figsize=(8,6),nrows=2,sharex=True,gridspec_kw={'height_ratios':[2,1]})
	fig.subplots_adjust(hspace=0.075)
	
	# formatting
	axs[0].plot(vals,mean_ngauss,'-o',color='b',label='Shared Area')
	# axs[0].plot(vals,mean_gauss,'--o',color='b')
	axs[0].plot(vals,mean_cc,'-s',color='r',label='Cross Corr')
	# axs[0].plot(vals,mean_ccgauss,'--s',color='r')
	axs[1].plot(vals,std_ngauss/mean_ngauss,'-o',color='b')
	# axs[1].plot(vals,std_gauss/mean_gauss,'--o',color='b')
	axs[1].plot(vals,std_cc/mean_cc,'-s',color='r')
	# axs[1].plot(vals,std_ccgauss/mean_ccgauss,'--s',color='r')

	# axs[0].errorbar(vals,mean_ngauss,yerr=std_ngauss,fmt='o',color='b')
	# axs[0].errorbar(vals,mean_gauss,yerr=std_gauss,fmt='s',color='g')
	
	plt.gca().ticklabel_format(useOffset=False)
	axs[0].set_xlim(0,vals[-1])
	axs[1].set_ylim(0,1)
	axs[0].axhline(a2-a1,linestyle='--',color='darkgrey')
	axs[0].legend(loc='best')
	axs[0].set_ylabel('Offset, '+r'$\tau$',**tnrfont)
	axs[1].set_ylabel(r'$\sigma_\tau/\mu_\tau$',**tnrfont)
	axs[1].set_xlabel('Noise amp./Signal amp.',**tnrfont)
	plt.savefig('correlation/noisy-shared-area_MC.png')
	plt.show()
	return None

if __name__=='__main__':
	# package and list_new load
	from func_load import *
	shared_area_noise()
