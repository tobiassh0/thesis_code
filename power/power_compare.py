import numpy as np
import matplotlib.pyplot as plt 
import my_constants as const
from list_new import *


###################
mD = const.me_to_mD
mT = const.me_to_mT
###################

# could make a class variable to handle this calculation, otherwise I'm just typing it out multiple times

#sim_lst = ['traceT_0_00','traceT_D_99_T_01','traceT_D_89_T_11','traceT_0_50']#'cold_JET26148','traceT_0_11']#,L_8rL','L_10rL','L_12rL','L_13rL','L_20rL','L_22rL','L_27rL']
sim_lst = ['traceT_0_50','traceT_D_89_T_11','traceT_D_99_T_01','traceT_0_00']
#'traceT_D_70_T_30','traceT_D_82_T_18','traceT_D_95_T_05',

#height = 5 ; width = 10
height = 3 ; width = 10

fig, axs = plt.subplots(nrows=1,figsize=(width,height),sharex=True)
fig.subplots_adjust(hspace=0.1)
labels = [r'$50\%$' +'  '+r'$T_3$',r'$11\%$'+'  '+r'$T_3$',r'$1\%$'+'  '+r'$T_3$',r'$0\%$' +'  '+r'$T_3$']#,'cold']
#,r'$30\%$' +'  '+r'$T_3$',r'$18\%$' +'  '+r'$T_3$',r'$5\%$' +'  '+r'$T_3$',
wmin = 10 ; wmax = 25

D_cyc = np.arange(0,int(wmax+1),1) ; mD_mT = const.me_to_mD/const.me_to_mT
T_cyc = D_cyc*mD_mT
for j in range(len(D_cyc)):
	axs.axvline(D_cyc[j] ,linestyle='--',color='darkgrey')

powArr = []
mu = np.array([50,30,18,11,5,1,0])
peak_freq = np.zeros(len(sim_lst))
for i in range(len(sim_lst)):
	sim_loc = getSimulation('/storage/space2/phrmsf/'+sim_lst[i])
	wnorm = getCyclotronFreq(sdfread(0),'Deuterons')
	times = read_pkl('times')
	power = 10**read_pkl('log10_power')
	omegas = read_pkl('omegas_power')
	if sim_lst[i] in ['traceT_0_00']:
		dw = wnorm*(omegas[-1]-omegas[0])/len(omegas)
	else:
		dw = (omegas[-1]-omegas[0])/len(omegas)
		omegas = omegas/wnorm
	psd = power/dw
	peak_freq[i] = omegas[np.argmax(psd)]
	axs.plot(omegas,psd,label=labels[i])
	powArr.append(psd)		

#p = np.polyfit(mu,peak_freq,deg=1) ; MU = np.arange(0,50,0.1)
#plt.plot(mu,peak_freq,'ko',MU,MU*p[0]+p[1],'--r')
#plt.xlabel(r'$\xi_T$',fontsize=18) ; plt.ylabel(r'$\omega_{max}/\Omega_D$',fontsize=18)
#plt.show() ; fig.savefig('/storage/space2/phrmsf/paper/mu_omega_max.png',bbox_inches='tight')

powArr = np.array(powArr)
## normal
axs.set_ylim(min(powArr[0,10:])/2,5*max(powArr[0,10:]))
#axs.legend(loc='best')
axs.set_yscale('log')
axs.set_xlim(wmin,wmax)
axs.set_ylabel(r'PSD',**tnrfont)
axs.set_xlabel(r'$\omega/\Omega_D$',fontsize=20)
fig.savefig('/storage/space2/phrmsf/paper/Power_log_zoom.png',bbox_inches='tight')
#fig.savefig('/storage/space2/phrmsf/paper/Power_log_zoom.eps',bbox_inches='tight')

## EPS 4-page
#axs.set_ylim(min(powArr[0,10:])/2,5*max(powArr[0,10:]))
#axs.legend(loc='best',frameon=False)
#axs.set_yscale('log')
#axs.set_xlim(wmin,wmax)
#axs.set_ylabel(r'Power Spectral Density',**tnrfont)
#axs.set_xlabel(r'Frequency  '+r'$[\Omega_D]$',**tnrfont)
#fig.savefig('/storage/space2/phrmsf/paper/EPS-4page/PSD_log.png',bbox_inches='tight',dpi=75)


plt.show()




################################################################################




#import numpy as np
#import matplotlib.pyplot as plt 
#import my_constants as const
#from list_new import *


####################
#mD = const.me_to_mD
#mT = const.me_to_mT
####################

## could make a class variable to handle this calculation, otherwise I'm just typing it out multiple times

#sim_lst = ['traceT_0_00','traceT_0_01','traceT_0_11']#'cold_JET26148','traceT_0_11']#,L_8rL','L_10rL','L_12rL','L_13rL','L_20rL','L_22rL','L_27rL']
#fig, axs = plt.subplots(nrows=2,figsize=(12,7),sharex=True)
#fig.subplots_adjust(hspace=0.1)
#labels = [r'$0\%$' +'  '+r'$T_3$',r'$1\%$'+'  '+r'$T_3$',r'$11\%$'+'  '+r'$T_3$']

#powArr = []
#for i in range(len(sim_lst)):
#	sim_loc = getSimulation('/storage/space2/phrmsf/'+sim_lst[i])
#	wnorm = getCyclotronFreq(sdfread(0),'Deuterons')
#	times = read_pkl('times')
#	power = 10**read_pkl('log10_power')
#	omegas = read_pkl('omegas_power')
#	dw = wnorm*(omegas[-1]-omegas[0])/len(omegas)
#	psd = power/dw
#	axs[0].plot(omegas,psd)
#	powArr.append(psd)		

#wmin = 0 ; wmax = 22
#D_cyc = np.arange(0,35,1) ; mD_mT = const.me_to_mD/const.me_to_mT
#T_cyc = D_cyc*mD_mT
#powArr = np.array(powArr)
#psd_labels=[r'$PSD_{1\%}/PSD_{0\%}$',r'$PSD_{11\%}/PSD_{0\%}$',r'$PSD_{11\%}/PSD_{1\%}$']
#axs[0].legend(loc='lower right',labels=labels)
#axs[0].set_ylim(min(powArr[1,10:])/2,5*max(powArr[1,10:]))
#axs[1].plot(omegas,powArr[1,:]/powArr[0,:],color='darkgoldenrod') # 1%/0%
#axs[1].plot(omegas,powArr[2,:]/powArr[0,:],color='darkcyan') # 11%/0%
#axs[1].plot(omegas,powArr[2,:]/powArr[1,:],color='mediumturquoise') # 11%/1%
#axs[1].legend(loc='lower center',labels=psd_labels,ncol=3)
#for i in range(len(axs)): 
#	axs[i].set_yscale('log')
#	axs[i].set_xlim(wmin,wmax)
#	for j in range(len(D_cyc)):
#		axs[i].axvline(D_cyc[j] ,linestyle='--',color='k',alpha=0.5)
#		axs[i].axvline(T_cyc[j],linestyle='-.',color='r',alpha=0.5)
#		if D_cyc[j] <= wmax-1: 
#			if D_cyc[j] % 2 ==0:
#				axs[i].annotate(str(np.around(D_cyc[j],1)),xy=[D_cyc[j]/wmax+0.01,0.8],xycoords='axes fraction',color='k')
#		if T_cyc[j] <= wmax: 
#			if np.around((T_cyc[j] % 2),1) ==0:
#				axs[i].annotate(str(np.around(D_cyc[j],1)),xy=[T_cyc[j]/wmax+0.01,0.87],xycoords='axes fraction',color='r')

#axs[0].set_ylabel(r'PSD',fontsize=18)
#axs[1].set_xlabel(r'$\omega/\Omega_D$',fontsize=18)
#axs[1].set_ylabel(r'PSD ratios',fontsize=18)
#fig.savefig('/storage/space2/phrmsf/Power_compare_log_v2.png',bbox_inches='tight')
##plt.show()


