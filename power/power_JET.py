
from func_load import *


def signedCurv(fx,dx):
	### take the absolute for defined "curvature" or the square integral to get the total curvature of fx 
	## In
	#	fx : np array of points to calculate curvature ; dx : spacing in x-array points
	## Out
	#	kx : signed curvature, as defined here (https://en.wikipedia.org/wiki/Curvature#Graph_of_a_function)
	#	kxint : Integral of squared signed-curvature, higher the value --> more curvy / less smooth
	fxp = np.gradient(fx,dx)
	fxpp = np.gradient(fxp,dx)
	kx = abs(fxpp)/(1+fxp**2)**(1.5)
	kxint = np.sum((kx)*dx)
	return kx, kxint
	
###################
mD = const.me_to_mD
mT = const.me_to_mT
###################

sim_lst = ['traceT_0_00','traceT_D_99_T_01','traceT_D_89_T_11','traceT_0_50']#'cold_JET26148','traceT_0_11']#,L_8rL','L_10rL','L_12rL','L_13rL','L_20rL','L_22rL','L_27rL']
height = 7.5 ; width = 6.5
fig, axs = plt.subplots(nrows=len(sim_lst)+1,figsize=(width,height),sharex=True)
fig.subplots_adjust(hspace=0.1)
labels = [r'$0\%$' +'  '+r'$T_3$',r'$1\%$'+'  '+r'$T_3$',r'$11\%$'+'  '+r'$T_3$',r'$50\%$' +'  '+r'$T_3$']

JETdata = np.loadtxt('JET26148_ICE_POWER.txt',delimiter=',')
JETpower, JETfreqs = JETdata[:,1], JETdata[:,0] # 2 columns, N rows
wnorm = 2*const.PI*17E6
print('wnorm [MHz] :: {}'.format(wnorm/(2*const.PI*1e6)))
labelAll =[r'$\#26148$',r'$0\%$',r'$1\%$',r'$11\%$',r'$50\%$']
JETfreqs = 2*const.PI*JETfreqs*1e6/(wnorm) # convert MHz to wcD
maxJETfreqs = round(max(JETfreqs))
print('MAX FREQ (JET wcD):: ',maxJETfreqs)

## interpolate data
Ninterp = 10000
freqs = np.linspace(min(JETfreqs),max(JETfreqs),Ninterp)
JETpower = np.interp(freqs,JETfreqs,JETpower) # smooth out power so not as jagged
df = (freqs[-1]-freqs[0])/len(freqs)
KX, KXint = signedCurv(JETpower,df)
#ax2 = axs[0].twinx()
#ax2.plot(freqs,KX,color='b') ; ax2.set_ylim(0,max(abs(KX**2))+max(abs(KX**2))/10) ; ax2.set_ylabel(r'$k^2(\omega)$',fontsize=18,color='b')
dw = wnorm*(freqs[-1]-freqs[0])/len(freqs)
bs_JET = (1/len(JETpower))*np.sum((np.mean(JETpower)-JETpower)**2)
print('jet smoothness:',KXint)

psd = True
powArr = []
smoothArr = []
BS = []
BSS = []
smoothArr.append(KXint)

Dcyc = np.arange(0,maxJETfreqs+1,1) # MHz
for i in range(len(axs)):
	axs[i].annotate(labelAll[i],xy=(0.025,0.7),xycoords='axes fraction',fontsize=18,color='r')
	axs[i].set_xlim(0,maxJETfreqs)
	if i == 0: 
		axs[i].set_ylabel(r'$P_{ICE}$'+'  '+'[dB]',fontsize=18)
	else:
		if psd:
			axs[i].set_ylabel('PSD',fontsize=18)
		else:
			axs[i].set_ylabel(r'$\rho(\omega)$',fontsize=18) # probability density
	for h in Dcyc: 
		axs[i].axvline(h,linestyle='--',color='darkgrey')

## plot JET data
axs[0].plot(freqs,JETpower,color='k')
axs[0].set_xlim(0,maxJETfreqs)
JETpower = JETpower/(dw*np.sum(JETpower))

for i in range(len(sim_lst)):
#	ax2 = axs[i+1].twinx()
	sim_loc = getSimulation('/storage/space2/phrmsf/'+sim_lst[i])
	wnorm = getCyclotronFreq(sdfread(0),'Deuterons')
	wcD = getCyclotronFreq(sdfread(0),'Deuterons')
	times = read_pkl('times')
	power = 10**read_pkl('log10_power')
	omegas = read_pkl('omegas_power')
	if sim_lst[i] in ['traceT_0_50','traceT_D_89_T_11','traceT_D_99_T_01']: # un-normalised omegas
		dw = (omegas[-1]-omegas[0])/len(omegas)
		omegas = omegas/wnorm
	else:
		dw = wnorm*(omegas[-1]-omegas[0])/len(omegas)
	print('dw ::',dw)
#	power = np.interp(freqs,omegas,power) # interpolate
#	dw = wnorm*(freqs[-1]-freqs[0])/len(freqs)
	kx, kxint = signedCurv(power,dw)
	smoothArr.append(kxint)
	if psd: # psd
		new_power = power/dw
	else: # probability density
		new_power = power/(dw*np.sum(power))
		print('SUM ::',np.sum(dw*new_power))
#	thresh = omegas < maxJETfreqs
#	new_power = new_power[thresh]
#	axs[i+1].plot(omegas[thresh],(new_power),label=labels[i],color='k')  ; axs[i+1].set_yscale('log') ; axs[i+1].set_ylim(10**(np.log10(min(new_power))-0.5),10**(np.log10(max(new_power))+0.5)) # offset axis # *wcD/(2*const.PI*1e6)
	axs[i+1].plot(omegas,np.log10(new_power),label=labels[i],color='k')#  ; axs[i+1].set_yscale('log')
	axs[i+1].set_ylim(.5,3.)
#	ax2.plot(omegas,kx,color='b') ; ax2.set_ylim(0,max(abs(kx))+max(abs(kx))/10)
	## baseline
	Pmean = np.mean(np.log10(new_power))
	Pdiff = np.std(np.log10(new_power))
	Pstd = np.abs(np.std(new_power))
#	baseline  = Pdiff*np.cos(freqs*2*const.PI)+Pmean
#	axs[i+1].plot(freqs,abs(baseline),color='b')
#	kxb, kxbint = signedCurv(baseline, dw)
#	smoothArr.append(str('basek :: '+str(kxbint)))
#	ax2.set_ylabel(r'$k(\omega)$',fontsize=18,color='b')

print(labelAll)
print(smoothArr)
#print(BS)
#print(BSS)
axs[len(axs)-1].set_xlabel(r'$\omega/\Omega_D$',fontsize=18)
#plt.show()
fig.savefig('/storage/space2/phrmsf/paper/Power_ICE_compare_normWcD.png')
fig.savefig('/storage/space2/phrmsf/paper/Power_ICE_compare_normWcD.eps')


