
from func_load import *

os.chdir('/storage/space2/phrmsf/lowres_D_He3')
sims = np.sort([i for i in os.listdir() if 'p_90' in i])
labels = [i[2:4] for i in sims]
# rearrange so biggest to smallest T-concentration
print(sims)
wmax = 23
B0 = 3.7 # T
wcp = const.qe*B0/getMass('Protons') # rad/s
colors = plt.cm.rainbow(np.linspace(0,1,len(sims)))

fig,ax=plt.subplots(figsize=(13,4))

for c in range(len(sims)):
	simloc = getSimulation(sims[c])
	omegas = read_pkl('omegas_power')
	psd = (10**read_pkl('log10_power'))/((omegas[-1]-omegas[0])/len(omegas)) # power / dw
	
	# thresh below wmax + 1 (to remove non-ICE signal)
	thresh = omegas/wcp < wmax+1
	omegas = omegas[thresh]
	psd = psd[thresh]
	#plt.plot(omegas/wcp,np.log10(psd))
	#for i in np.arange(0,wmax+1,1):
	#	plt.axvline(i,linestyle='--',color='darkgrey')

	# find peaks
	peaks = extractPeaks(np.log10(psd),prominence=0.2) # hard-coded for now as peaks are messy (resolution) 
	#plt.scatter(omegas[peaks]/wcp,np.log10(psd[peaks]),color='k')
	#plt.show()
	ppsd = psd[peaks]
	pomegas = omegas[peaks]

	# empty arrays to plot line of abs(difference)
	sl = []
	sdw = []	

	# loop through harmonics
	l = np.arange(0,wmax,1)
	for i in range(len(l)):
		diffs_w = pomegas - l[i]*wcp # array of differences
		diffs_w = [i for i in diffs_w if (abs(i) < 0.5*wcp)] # less than half a harmonic away
		for dw in diffs_w:
			if dw != 0:
				sdw.append(dw)
				sl.append(l[i])
				#plt.scatter(l[i],abs(dw)/wcp,color=colors[c])
			else:
				continue
	ax.plot(sl,np.abs(sdw)/wcp,'o-',color=colors[c],label=labels[c])	
	os.chdir('..') # back to home

ax.legend(loc='best')
ax.set_ylim(-0.2,0.8)
ax.set_ylabel(r'$|\omega_{l} - l\Omega_\alpha|/\Omega_\alpha$',**tnrfont)
ax.set_xlabel(r'$l$',**tnrfont)
plt.show()
os.chdir('/storage/space2/phrmsf/lowres_D_He3')
fig.savefig('dw_l.png',bbox_inches='tight')
