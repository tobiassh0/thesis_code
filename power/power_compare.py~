
from func_load import *

os.chdir('/storage/space2/phrmsf/lowres_D_He3/')
sims = np.sort([i for i in os.listdir() if 'p_90' in i])
hlabels = np.array([int(i[2:4]) for i in sims])
power_compare(sims,labels=hlabels,normspecies='Deuterons')
sys.exit()

###############################################################################
###############################################################################

sim_lst = ['traceT_0_50','traceT_D_89_T_11','traceT_D_99_T_01','traceT_0_00']#,'cold_JET26148']
#'traceT_D_70_T_30','traceT_D_82_T_18','traceT_D_95_T_05',

height = 5 ; width = 10
#height = 3 ; width = 10

fig, axs = plt.subplots(nrows=1,figsize=(width,height),sharex=True)
fig.subplots_adjust(hspace=0.1)
labels = [r'$50\%$',r'$11\%$',r'$1\%$',r'$0\%$']#,'cold']
#,r'$30\%$' +'  '+r'$T_3$',r'$18\%$' +'  '+r'$T_3$',r'$5\%$' +'  '+r'$T_3$',
wmin = 0 ; wmax = 25

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
	dw = (omegas[-1]-omegas[0])/len(omegas)
	omegas = omegas/wnorm
	psd = power/dw
	peak_freq[i] = omegas[np.argmax(psd)]
	axs.plot(omegas,psd,label=labels[i])
	powArr.append(psd)

powArr = np.array(powArr)
## normal
axs.set_ylim(min(powArr[0,10:])/2,5*max(powArr[0,10:]))
axs.legend(loc='best')
axs.set_yscale('log')
axs.set_xlim(wmin,wmax)
axs.set_ylabel(r'PSD',**tnrfont)
axs.set_xlabel(r'$f/f_{cD}$',fontsize=20)
fig.savefig('/storage/space2/phrmsf/paper/remake/Power_log.png',bbox_inches='tight')
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
