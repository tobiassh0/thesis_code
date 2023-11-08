
from func_load import *

###################
mD = const.me_to_mD
mT = const.me_to_mT
###################

sim_lst = ['traceT_0_00','traceT_D_99_T_01','traceT_D_95_T_05','traceT_D_89_T_11','traceT_D_82_T_18','traceT_D_70_T_30','traceT_0_50']
height = 9 ; width = 7
fig, axs = plt.subplots(nrows=2,figsize=(width,height))#,sharex=True)
fig.subplots_adjust(hspace=0.3)
labels = [r'$0\%$' +'  '+r'$T_3$',r'$1\%$'+'  '+r'$T_3$',r'$5\%$'+'  '+r'$T_3$',r'$11\%$'+'  '+r'$T_3$',r'$18\%$'+'  '+r'$T_3$',r'$30\%$'+'  '+r'$T_3$',r'$50\%$'+'  '+r'$T_3$']
colors = list(reversed(['blue','pink','orange','g','k','r','cyan']))
tcolors = ['m','orange','g','royalblue']
freqoff = np.zeros(len(sim_lst))
zeroperc = False
Ninterp = 10000
tfreqoff = []
scolor = []
xi2 = [0,0.01,0.05,0.11,0.18,0.3,0.5]
CrossCorr = True ; PowerTrend = False ; SharedArea = False ; ArgminWoff = False
wcD = const.qe*2.1/(getMass('Deuterons'))

for i in range(len(sim_lst)):
	getSimulation('/storage/space2/phrmsf/'+sim_lst[i])
#	wcD = getCyclotronFreq(sdfread(0),'Deuterons')
	freq = read_pkl('omegas_power')
	power = 10**read_pkl('log10_power')
	tdw = (freq[-1]-freq[0])/len(freq)
	if sim_lst[i] in ['traceT_0_00']:
		freq = wcD*freq
		zeroperc = True
	else:
		zeroperc = False
#	print(sim_lst[i],max(freq))
	freqs = np.linspace(min(freq),max(freq),Ninterp) # interpolate so same length
	dw = (freqs[-1]-freqs[0])/len(freqs)
	psd = np.interp(freqs,freq,power)/(dw)#*np.sum(power0))
	if zeroperc:
		power0 = power
		psd0 = psd
		fft_power0 = np.fft.fft(power0)

#	## Cross Correlation
	if zeroperc: # auto
		corrArr = signal.correlate(psd0,psd0,mode='full')
		corrArrMax = np.max(corrArr)
		corrArr /= corrArrMax
	else:
		corrArr = signal.correlate(psd,psd0,mode='full')	
		corrArr /= corrArrMax
	tfreqs = freqs.copy()
	tfreqs = np.append(tfreqs,-tfreqs[1:])
	tfreqs = np.sort(tfreqs)
	axs[0].plot(tfreqs/wcD,corrArr,label=labels[i],color=colors[i])
	arg = np.argmax(corrArr)
	freqoff[i] = tfreqs[arg]
	axs[0].scatter(tfreqs[arg]/wcD,corrArr[arg],marker='o',edgecolors='k',c=colors[i],s=50)
	fcc = tfreqs[arg]

	## Trend in maximum power spectra
#	axs[0].plot(freqs/wcD,psd)
	arg = np.argmax(psd)
#	axs[0].scatter(freqs[arg]/wcD,psd[arg],c=colors[i],edgecolors='k',s=50,marker='o')
	freqoff[i] = freqs[arg]
	fmp = freqs[arg]
	if zeroperc:
		fmp0 = fmp
		fmp = 0
	else:
		fmp = fmp - fmp0
#	axs[1].scatter(xi2[i],fmp/wcD,color=tcolors[1])
#	
	## Phase correlation
	import phasecorrelation as pc
	shift =  pc.phaseCorrelation(power,fft_power0,tdw,wcD)
	fpc = -shift

	## Shared area between curves & Argmin difference
	tarea = np.zeros(len(psd))
	diff = np.zeros(len(psd))
	woff = np.zeros(len(psd))
	for r in range(len(psd)):
		woff[r] = (freqs[0] + (r-len(psd))*dw)
		tpsd = np.roll(psd,len(psd)-r)
		argsim_0 = tpsd >= psd0
		arg0_sim = argsim_0 == False
		tarea[r] = np.sum((tpsd*arg0_sim + psd0*argsim_0)*dw)
		diff[r] = np.sum(np.abs(tpsd-psd0))
	sarg = np.argmax(tarea)
	axs[0].plot(woff/wcD,tarea)
	axs[0].scatter(woff[sarg]/wcD,tarea[sarg],edgecolors='k',c=colors[i],marker='o')
	darg = np.argmin(diff)
	freqoff[i] = freqs[sarg]
	freqoff[i] = freqs[darg]
	if freqoff[i] > 17.5*wcD:
		freqoff[i] -= 35*wcD
#	axs[0].plot(woff/wcD,diff)
	fsa = freqs[sarg]
	fam = freqs[darg]
	if fsa > 17.5*wcD:
		fsa -= 35*wcD
	if fam > 17.5*wcD:
		fam -= 35*wcD
	axs[1].scatter(xi2[i],fsa/wcD,color=tcolors[2])
	axs[1].scatter(xi2[i],fam/wcD,color=tcolors[3])
	print(wcD)
	print(sim_lst[i],fcc,fmp,fsa,fpc)#,fam
	tfreqoff.append([fcc,fmp,fsa,fpc])#,fam]
	
##########################################################
##########################################################

if CrossCorr:
	axs[0].axvline(0,linestyle='--',color='darkgrey')
	axs[0].set_xlim(-6,2)
	axs[0].set_ylim(-0.01,2)
	axs[0].set_xlabel(r'$\omega_{off}/\Omega_D$',fontsize=18)
	axs[0].set_ylabel(r'Normalised cross-correlation',fontsize=18)
	axs[0].legend(loc='best')
	figname = 'power_crosscor'

if PowerTrend:
	freqoff = freqoff - freqoff[0]
	figname = 'power_trend'	
	axs[0].set_xlim(10,22)
	axs[0].set_ylim(10**2,10**9)
	axs[0].set_yscale('symlog')
	axs[0].set_ylabel(r'PSD',fontsize=18)
	axs[0].set_xlabel(r'$\omega/\Omega_D$',fontsize=18)

if SharedArea:
	figname = 'power_sharedarea'
	axs[0].set_xlim(-35,0)
	axs[0].set_yscale('symlog')
	axs[0].set_ylabel(r'Shared area',fontsize=18)
	axs[0].set_xlabel(r'$\omega_{off}/\Omega_D$',fontsize=18)

if ArgminWoff:
	figname = 'power_argmin'
	axs[0].set_xlim(-35,0)
	axs[0].set_ylabel(r'Minimum difference',fontsize=18)
	axs[0].set_xlabel(r'$\omega_{off}/\Omega_D$',fontsize=18)

#xi2 = [0,0.01,0.05,0.11,0.18,0.3,0.5]
#p = np.polyfit(xi2,freqoff,deg=1) ; XI2 = np.arange(0,0.50,0.001)
#print('PARAMETERS :: ',p)
#axs[1].plot(xi2,freqoff/wcD,'ko',XI2,(XI2*p[0]+p[1])/wcD,'--r')
#axs[1].set_ylabel(r'$\omega_{off}/\Omega_D$',fontsize=18)
#axs[1].set_xlabel(r'$\xi_T$',fontsize=18)
#axs[1].set_xlim(0,0.5)
#fig.savefig('/storage/space2/phrmsf/paper/power-comparison/'+figname+'.png',bbox_inches='tight')
#fig.savefig('/storage/space2/phrmsf/paper/power-comparison/'+figname+'.eps',bbox_inches='tight')


### compare growths and intercepts
## cross cor
#cc = [-4.69575316e+08, -2.20531068e+07]
## trend in max 
#tm = [-4.28126024e+08,  7.90327042e+06]
## shared area
#sa = [-5.18487095e+08, -6.50436519e+06] 
## argmin
#am = [-5.22885535e+08, -1.01009829e+07] 

#grads = np.array([cc[0],tm[0],sa[0],am[0]])/wcD
#intercepts = np.array([cc[1],tm[1],sa[1],am[1]])/wcD
#grad = np.mean(grads)
#intercept = np.mean(intercepts)
#stdgrad = np.std(grads)
#stdintercept = np.std(intercepts)

#print(str(grad)+' +- '+str(stdgrad))
#print(str(intercept)+' +- '+str(stdintercept))

####################################################

tfreqoff=np.array(tfreqoff,dtype='float')
fcc = tfreqoff[:,0]
fmp = tfreqoff[:,1]
fsa = tfreqoff[:,2]
fpc = tfreqoff[:,3]
#fam = tfreqoff[:,2]
xi2 = [0,0.01,0.05,0.11,0.18,0.3,0.5]

### Fitting deg=1 poly each and total
fig,ax = plt.subplots(figsize=(7,4))
XI2 = np.arange(0,0.51,0.001)
## each fit
fitcc = np.polyfit(xi2,fcc,deg=1)
fitmp = np.polyfit(xi2,fmp,deg=1)
fitsa = np.polyfit(xi2,fsa,deg=1)
fitpc = np.polyfit(xi2,fpc,deg=1)
## cross-corr
ax.plot(XI2,(XI2*fitcc[0]+fitcc[1])/wcD,color=tcolors[0])
## max power
ax.plot(XI2,(XI2*fitmp[0]+fitmp[1])/wcD,color=tcolors[1])
##shared area
ax.plot(XI2,(XI2*fitsa[0]+fitsa[1])/wcD,color=tcolors[2])
## phase corr
ax.plot(XI2,(XI2*fitpc[0]+fitpc[1])/wcD,color=tcolors[3])
### total
meangrad = np.mean([fitcc[0],fitmp[0],fitsa[0],fitpc[0]])#,fitam[0]
stdgrad  = np.std([fitcc[0],fitmp[0],fitsa[0],fitpc[0]])#,fitam[0]
meanintercept = np.mean([fitcc[1],fitmp[1],fitsa[1],fitpc[1]])#,fitam[1]
stdintercept = np.std([fitcc[1],fitmp[1],fitsa[1],fitpc[1]])#,fitam[1]
ax.plot(XI2,(meangrad*XI2+meanintercept)/wcD,'k--')
ax.fill_between(XI2,((meangrad+stdgrad)*XI2+(meanintercept+stdintercept))/wcD,((meangrad-stdgrad)*XI2+(meanintercept-stdintercept))/wcD,color='lightcoral') # biggest error range

print(meangrad/wcD, ' +- ',stdgrad/wcD)
print(meanintercept/wcD, ' +- ',stdintercept/wcD)

## Scatter plot each point
tfreqoff = np.reshape(tfreqoff,-1)
mean_methods = np.mean(tfreqoff.reshape(-1,4),axis=1)
std_methods = np.std(tfreqoff.reshape(-1,4),axis=1)
#ax.scatter(xi2,mean_methods/wcD,color='k')
ax.errorbar(xi2,mean_methods/wcD,yerr=std_methods/wcD,fmt='s',color='k')
methods = ['Cross-correlation','Peak trend','Shared area','Phase-correlation']
ax.legend(methods,loc='best',frameon=True,shadow=True)
ax.set_ylabel(r'$\omega_{off}/\Omega_D$',fontsize=18)
ax.set_xlabel(r'$\xi_T$',fontsize=18)
ax.set_xlim(0,0.51)
#ax.set_ylim(-2.8,0.2)
#fig.savefig('/storage/space2/phrmsf/paper/power-comparison/power_trend.png',bbox_inches='tight')
fig.savefig('/storage/space2/phrmsf/paper/power-comparison/power_trend.tiff',bbox_inches='tight')
#fig.savefig('/storage/space2/phrmsf/paper/power-comparison/power_trend.eps',bbox_inches='tight')

plt.show() ; sys.exit()






