
from func_load import *
from sklearn.neighbors import KernelDensity



def loadJETdata():
	print('loading JET data')
	JETdata = np.loadtxt('JET26148_ICE_POWER.txt',delimiter=',')
	JETpower, JETfreqs = JETdata[:,1], JETdata[:,0] # 2 columns, N rows
	wnorm = 2*const.PI*17E6
	JETfreqs = 2*const.PI*JETfreqs*1e6/(wnorm) # convert MHz to wcD
	maxJETfreqs = round(max(JETfreqs))
	return maxJETfreqs, JETfreqs, JETpower

def loadSimData(wlim):
	power = read_pkl('log10_power')
	wcD = getCyclotronFreq(sdfread(0),'Deuterons')
	omegas = read_pkl('omegas_power')/wcD # normalised
#	## limit to a certain frequency
#	maxFreq = wlim #round(max(omegas))
#	thresh = omegas < maxFreq
#	omegas = omegas[thresh] ; power = power[thresh]

	return omegas, power

def funcparabola(X,c0,c1,c2):
	return c2+c1*(c0-X)**2

def funcexp(X,c0,c1,c2):
	return np.exp(c0*(X-c1))-c0*(X-c1)-c2

def funcnegexp(X,c0,c1,c2):
	return c0*(X-c1) + np.exp(-c0*(X-c1)) + c2

def funcquart(X,c0,c1,c2,c3,c4):
	return c0*X**4 + c1*X**3 + c2*X**2 + c3*X**1 + c4

def extractPeaks(power, Nperwcd=1):
	return signal.find_peaks(power,distance=Nperwcd)[0] # tune till Nperwcd encapsulates all peaks (visually)

def getPowerRatios(power,peaks):
	ratios = np.zeros(len(peaks))
	for i in range(1,len(peaks)):
		ratios[i] = power[peaks[i]]/power[peaks[i-1]]
	return ratios

def calculateTauSquared(x_array,y_array,matrix_norm,xbins,ybins,modelProb=0.7): 

	# Calculates tau squared given input data and model. 
	# The floor value is defined when making the model matrix. The floor is 
	# assumed to be a uniform probability across the period mass plane. 
	# in: 
	#   x_array: array of observed x data: mass in Msun 
	#   y_array: array of observed y data: either log or linear rotation period 
	#   matrix_norm: array of normalised model probability, either log or linear ## KDE of y-axis 
	#   ymax: maximum extent of rotation period, in either linear or log space 
	#   ymin: minimum extent of rotation period, in either linear or log space 
	#   modelProb: defines a model (and floor) probability as a fraction = 0.7 
	# out: 
	#   tau_2: the total tau 2 value of the model and the data 
	# To calculate tau squared, loop over data 

	tau_2 = 0. 
	xmin = 0; xmax = max(x_array)
	ymin = -3; ymax = 3 
	# makes NxM matrix of tau-squared 
	tau_2_arr = []
	for i in range(len(x_array)): 
		# Find which bin the data goes into 
		x_bin = int((x_array[i] - xmin)*xbins/(xmax - xmin)) ## N = no. bins in x 
		y_bin = int((y_array[i] - ymin)*ybins/(ymax - ymin)) ## M = no. bins in y 
		# Find the model density of that bin
		density_cell = matrix_norm[x_bin][ybins-1-y_bin]
		# calculate the (model + background) density rhoPrime 
		rhoPrime = (1.- modelProb)/(xbins*ybins*(xmax-xmin)*(ymax-ymin)) + modelProb*density_cell # TODO; how to calculate area
		# =========Calculate tau squared ============= 
		tau_2 += -2.*np.log(rhoPrime) 
		tau_2_arr.append(-2.*np.log(rhoPrime))
	return tau_2, tau_2_arr 


#fig, ax = plt.subplots(figsize=(10,5.5))#,subplot_kw=dict(projection='polar'))
fig, ax = plt.subplots(figsize=(8,5))#,subplot_kw=dict(projection='polar'))
## load JET data
maxJETfreqs, JETfreqs, JETpower = loadJETdata()
JETpeaks = extractPeaks(JETpower,1)
ratiosJET = getPowerRatios(JETpower, JETpeaks)

## calculate KDEs of JET data
Nval = 2000
xbins = maxJETfreqs # same number excluding 0th bin, centred on the harmonic
KDEs = np.zeros((xbins,Nval))
y_d = np.linspace(-3,3,Nval)
for i in range(xbins):
	# instantiate and fit the KDE model
	kde = KernelDensity(bandwidth=0.3, kernel='gaussian') # fixed bandwidth (taken from previous work on Stellar clusters)
	# find values of ratiosJET that fit within a given bin
	bin_range = [i-0.5, i+0.5]
	val = ratiosJET[(JETfreqs[JETpeaks][:]>bin_range[0]) & (JETfreqs[JETpeaks][:]<bin_range[1])] # normalised freqs
	pp = (JETfreqs[JETpeaks][:]>bin_range[0]) & ((JETfreqs[JETpeaks][:]<bin_range[1]))
	J = (JETfreqs[JETpeaks])
#	print((J[pp]))
	val = [[i] for i in val] # make N dimensional with N = len(val)
	kde.fit(val)
	# score_samples returns the log of the probability density
	logprob = kde.score_samples(y_d[:, None])
#	print(np.sum(np.exp(logprob)*(x_d[-1]-x_d[0])/len(x_d))) # should sum to 1
#	plt.fill_between(x_d, np.exp(logprob), alpha=0.5)
#	plt.plot(ratiosJET, np.full_like(ratiosJET, -0.01), '|k', markeredgewidth=1)
#	plt.ylim(-0.02, 0.22)
	print(np.sum(np.exp(logprob)*(y_d[-1]-y_d[0])/len(y_d))) # makes sure prob along given freq bin sums to 1
	KDEs[i,:] = np.exp(logprob)/xbins

print('total KDE sum here:',np.sum(KDEs*(y_d[-1]-y_d[0])/len(y_d))) # makes sure prob along given freq bin sums to 1

## load sim data
sims = ['traceT_0_00','traceT_D_99_T_01','traceT_D_95_T_05','traceT_D_89_T_11','traceT_D_82_T_18','traceT_D_70_T_30','traceT_0_50']
c = 0
markers=['x','o','^','d','>','2','*']
labels = [r'$0\%$',r'$1\%$',r'$5\%$',r'$11\%$',r'$18\%$',r'$30\%$',r'$50\%$']
tau2_Npeak=[]
tau2_arrs=[]
tau2bar=[]
skip=0
for sim in sims:
	sim_loc = getSimulation('/storage/space2/phrmsf/'+sim)
	omegas, power = loadSimData(wlim=maxJETfreqs)
	wcD = getCyclotronFreq(sdfread(0),'Deuterons')
	if sim in ['traceT_0_00','traceT_0_11','traceT_0_01']:
		omegas = omegas*wcD
	## limit to a certain frequency
	maxFreq = maxJETfreqs #round(max(omegas))
	thresh = omegas < maxFreq
	omegas = omegas[thresh] ; power = power[thresh]
	## extract peaks
	peaks = extractPeaks(power)
#	plt.scatter(omegas[peaks],power[peaks])
#	plt.plot(omegas,power,label=labels[c]) ; plt.show()
	## ratios
	ratios = getPowerRatios(power, peaks)
	tau2, tau2arr = calculateTauSquared(omegas[peaks],ratios,KDEs,xbins-1,len(y_d))
	tau2_arrs.append(np.array(tau2arr))#tau2/len(ratios))
	tau2 = np.sum(tau2arr[skip:])
	tau2_Npeak.append(tau2/(len(tau2arr[skip:])))
	tau2bar.append(np.mean(tau2arr[skip:])/len(tau2arr[skip:]))
	## plot
#	plt.scatter(omegas[peaks],ratios,color='k',marker=markers[c],label=labels[c])
	print(np.around(omegas[peaks],1))
	c+=1

colors = ['cyan','r','k','g','orange','orchid','b']
#tau2_arrs = np.array(tau2_arrs,dtype='object').reshape((len(sims),-1))
#print(tau2_arrs.shape,tau2_arrs[0][1])

### show tau2 heatmap
#plt.xlabel(r'$\omega/\Omega_D$',fontsize=20)
#plt.ylabel(r'$PPR$',fontsize=20)
#im = plt.imshow(KDEs.T,cmap='jet',aspect='auto',origin='lower',interpolation='none',extent=[-0.5,10.5,min(y_d),max(y_d)])
#cbar = plt.colorbar(im).set_label(label=r'$\rho(\omega,r)$',size=20)
#plt.ylim(0.,3) ; plt.xlim(0,10.5)
#plt.legend(loc='best')
#plt.show()
#fig.savefig('/storage/space2/phrmsf/dump/tau_squared.eps',bbox_inches='tight')

### show tau2_Npeak per Mu
#print(tau2_Npeak)
#print(len(tau2bar))
#Mu = [0,1,5,11,18,30,50]
##plt.scatter(Mu,tau2_Npeak,color='k')
##plt.ylabel(r'$\tau^2/N_{peaks}$',fontsize=20)
#print(tau2bar)
#plt.ylabel(r'$\bar{\tau_i^2} /N_{peaks}$',fontsize=20)
#plt.xlabel(r'$\mu$'+'  '+'[%]',fontsize=20)
##popt, pcov = curve_fit(funcparabola,Mu,tau2_Npeak)
##popt, pcov = curve_fit(funcnegexp,Mu,tau2bar,bounds=([-1,5,-5],[0,20,10]),maxfev=10000)
#popt, pcov = curve_fit(funcquart,Mu,tau2bar,maxfev=10000)#,bounds=([-1,5,-5],[0,20,10]),maxfev=10000)
#perr = np.sqrt(np.diag(pcov))
#print(popt,perr)
#muarr = np.arange(-10,60,0.1)
##ytau = funcnegexp(muarr,popt[0],popt[1],popt[2])
#ytau = funcquart(muarr,popt[0],popt[1],popt[2],popt[3],popt[4])
#plt.plot(muarr,ytau,color='k',linestyle='--')
#plt.xlim(-10,60) ; plt.ylim(1.8,2.6)
#plt.fill_between(muarr,ytau-perr[-1],ytau+perr[-1],color='lightcoral')#,alpha=0.5)
#plt.scatter(Mu,tau2bar,color='k')
#fig.savefig('/storage/space2/phrmsf/paper/tau2bar_Mu_fit_quart.eps'.format(skip),bbox_inches='tight')
#fig.savefig('/storage/space2/phrmsf/paper/tau2bar_Mu_fit_quart.png'.format(skip),bbox_inches='tight')
#plt.show()

#### show tau2 contribution per peak
#tau2_arrs = np.array(tau2_arrs,dtype='object')
#for i in range(tau2_arrs.shape[0]):
#	index = np.arange(len(tau2_arrs[i]))
#	plt.scatter(index,tau2_arrs[i]/len(index),s=50,marker='x',c=colors[i],label=labels[i])
#	print(labels[i],tau2_arrs[i]/len(index))
#plt.legend(loc='center left',bbox_to_anchor=(1.,0.5))
#plt.ylim(1.95,2.6)
#plt.xlim(0,14.5)
#plt.ylabel(r'$\tau^2_i/N_{peaks}$',fontsize=20)
#plt.xlabel(r'$PPR(i)$',fontsize=20)
#fig.savefig('/storage/space2/phrmsf/paper/tau2_Npeak.eps',bbox_inches='tight')
#fig.savefig('/storage/space2/phrmsf/paper/tau2_Npeak.png',bbox_inches='tight')
#plt.show()


### plot parabolas per set of peaks
tau2_arrs = np.array(tau2_arrs,dtype='object')
# convert tau2 array to 2d array with None for missing peaks 
maxpeaks = 0
for i in range(tau2_arrs.shape[0]):
	tau2_arrs[i] = tau2_arrs[i]/len(tau2_arrs[i]) # tau2 per-peak
	if len(tau2_arrs[i]) > maxpeaks:
		maxpeaks = len(tau2_arrs[i])
	else:
		continue
print(maxpeaks)
ttau2_arrs = np.empty((len(sims),maxpeaks))
ttau2_arrs.fill(np.nan)
for i in range(tau2_arrs.shape[0]):
	j = 0
	while j < maxpeaks+1:
		try:
			ttau2_arrs[i,j] = tau2_arrs[i][j]
			j+=1
		except:
			j = maxpeaks+2
#print(ttau2_arrs)

fig = plt.figure(figsize=(10,10)) ; ax = plt.axes(projection='3d')
MuMin=[]
perr_arr=[]
zmin = 1.6 ; zmax = 2.6
# read each peak (column) and find parabola
for j in range(2,ttau2_arrs.shape[1]):
	tau2_parabola=np.zeros(len(sims))
	calculate=True
	for i in range(ttau2_arrs.shape[0]):
		if not np.isnan(ttau2_arrs[i,j]):
			tau2_parabola[i] = ttau2_arrs[i,j]
		else:
			calculate=False
		### 
	if calculate:
		mu = [0,1,5,11,18,30,50]
		popt,pcov=curve_fit(funcnegexp,mu,tau2_parabola,maxfev=200000)
		perr = np.sqrt(np.diag(pcov))
		MuMin.append(popt[1])
		perr_arr.append(perr[1])
#		print(popt,perr)
		## plot in 3d
		print(np.mean(tau2_parabola))
		PEAK = np.ones(len(tau2_parabola))*j
		ax.scatter(PEAK,mu,tau2_parabola,color=colors)
		Mu = np.arange(0,50,1)
		PEAKS = np.ones(len(Mu))*j
#		ax.plot(PEAKS,Mu,funcnegexp(Mu,popt[0],popt[1],popt[2]),color='k',linestyle='-')
#		ax.plot([j,j],[popt[0],popt[0]],[0,popt[2]],color='k',linestyle='--')
#		ax.plot([j,j],[popt[1],popt[1]],[zmin,funcnegexp(popt[1],popt[0],popt[1],popt[2])],color='k',linestyle='--')
ax.set_xlabel(r'$PPR(i)$',fontsize=20)
ax.set_ylabel(r'$\mu$'+'  '+r'$[\%]$',fontsize=20)
ax.set_zlabel(r'$\tau_i^2/N_{peaks}$',fontsize=20)
ax.set_zlim(zmin,zmax)
ax.set_ylim(-0.1,55)
ax.view_init(elev=15, azim=-145)
plt.show()
#fig.savefig('/storage/space2/phrmsf/paper/tau2_exp.png',bbox_inches='tight')
#fig.savefig('/storage/space2/phrmsf/paper/tau2_exp.eps',bbox_inches='tight')

perr_arr = np.array(perr_arr)
MuMin = np.array(MuMin)#,dtype='object')
d_2 = np.sum((MuMin**2)/perr_arr**2)/np.sum(1/(perr_arr**2))
d2 = (np.sum(MuMin/perr_arr**2)/np.sum(1/(perr_arr**2)))**2
Neff = (np.sum(1/perr_arr)**2)/(np.sum(1/perr_arr**2))
stderr = np.sqrt((d_2-d2)/Neff)

stderr = (np.sqrt(np.sum(perr_arr**2)))/len(perr_arr)
#print(np.mean(MuMin),' +- ',stderr)
  

