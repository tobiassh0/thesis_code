

def shared_area_noise():
	## shared area example
	
	# x-dimension array (arb)
	t = np.linspace(0,10,1000)
	f = 2*const.PI
	
	# spline fitting
	knot_numbers = 20
	
	# noise amplitude values to loop through
	vals = [0.,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0,3.25,3.5]
	
	# number of runs to do at each noise amplitude
	N=1e4
	peak_area=np.zeros((len(vals),int(N)))
	peak_spline=np.zeros((len(vals),int(N)))
	
	# signal parameters
	a1 = 1. ; b1 = 1
	a2 = 6. ; b2 = b1
	
	# Monte Carlo method
	for j in range(len(vals)):
		noise_floor = vals[j]*b1
		pa = np.zeros(int(N)) ; ps = np.zeros(int(N))
		for i in range(int(N)):
			# noise and signals	
			noise1 = np.random.normal(0,noise_floor,len(t)) + 5 # creates random noise with each iteration
			noise2 = np.random.normal(0,noise_floor,len(t)) + 5
			sig1 = np.exp(-b1*(t-a1)**2)+noise1
			sig2 = np.exp(-b2*(t-a2)**2)+noise2
	
			# shared area between signals
			areas = shared_area(sig1,sig2) # from list_new
	
			# fitting spline
			t_new = np.linspace(0, 1, knot_numbers+2)[1:-1]
			q_knots = np.quantile(t, t_new)
			tf,c,k = interpolate.splrep(t, areas, t=q_knots, s=1)
			yfit = interpolate.BSpline(tf,c,k)(t)
	
			# summing total
			pa[i]=t[np.argmax(areas)]
			ps[i]=t[np.argmax(yfit)]
		
		# append all values in 2d array
		peak_area[j,] = (pa)
		peak_spline[j,] = (ps)
	
	# mean and std of each run
	mean_area=np.mean(peak_area,axis=1)
	mean_spline=np.mean(peak_spline,axis=1)
	std_area=np.std(peak_area,axis=1)
	std_spline=np.std(peak_spline,axis=1)
	
	# setup plot
	fig,axs=plt.subplots(figsize=(10,10),nrows=2,sharex=True,gridspec_kw={'height_ratios':[3,1]})
	fig.subplots_adjust(hspace=0.075)
	
	# formatting
	axs[0].plot(vals,mean_area,'-o',vals,mean_spline,'--s')
	axs[1].plot(vals,std_area/mean_area,color='b')
	axs[1].plot(vals,std_spline/mean_spline,color='g')
	#axs[0].errorbar(vals,mean_area,yerr=std_area,fmt='o',color='b')
	#axs[0].errorbar(vals,mean_spline,yerr=std_spline,fmt='s',color='g')
	
	plt.gca().ticklabel_format(useOffset=False)
	axs[0].set_ylabel('Offset, '+r'$\tau$',**tnrfont)
	axs[1].set_ylabel(r'$\sigma_\tau/\mu_\tau$',**tnrfont)
	axs[1].set_xlabel('Noise amp./Signal amp.',**tnrfont)
	#plt.savefig('correlation/noisy-shared-area_MC_errorbars-test1.png')
	plt.show()
	return None

if __name__=='__main__':
	shared_area_noise()
	# package and list_new load
	from func_load import *
