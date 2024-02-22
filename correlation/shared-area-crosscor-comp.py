
from func_load import *
from matplotlib.gridspec import GridSpec
import correlation.phasecorrelation as pc

# total figure
fig = plt.figure(figsize=(8,8))
# base fig (raw data w/ noise)
gs1 = GridSpec(5,12,bottom=0.05,top=0.95,left=0.05,right=0.48,wspace=0.125,hspace=0.0)
ax00 = fig.add_subplot(gs1[0,:])
ax10 = fig.add_subplot(gs1[1,:],sharex=ax00)
ax20 = fig.add_subplot(gs1[2,:],sharex=ax00)
ax30 = fig.add_subplot(gs1[3,:],sharex=ax00)
ax40 = fig.add_subplot(gs1[4,:],sharex=ax00)
# each method sub-panel
gs2 = GridSpec(5,12,bottom=0.05,top=0.95,left=0.55,right=0.98,wspace=0.05,hspace=0.05)
ax01 = fig.add_subplot(gs2[0, :6])
ax02 = fig.add_subplot(gs2[0, 6:],sharey=ax01)
ax11 = fig.add_subplot(gs2[1, :6],sharex=ax01,sharey=ax01)
ax12 = fig.add_subplot(gs2[1, 6:],sharex=ax02,sharey=ax01)
ax21 = fig.add_subplot(gs2[2, :6],sharex=ax01,sharey=ax01)
ax22 = fig.add_subplot(gs2[2, 6:],sharex=ax02,sharey=ax01)
ax31 = fig.add_subplot(gs2[3, :6],sharex=ax01,sharey=ax01)
ax32 = fig.add_subplot(gs2[3, 6:],sharex=ax02,sharey=ax01)
ax41 = fig.add_subplot(gs2[4, :6],sharex=ax01,sharey=ax01)
ax42 = fig.add_subplot(gs2[4, 6:],sharex=ax02,sharey=ax01)

# combine into rows
axrow0 = [ax00,ax01,ax02]
axrow1 = [ax10,ax11,ax12]
axrow2 = [ax20,ax21,ax22]
axrow3 = [ax30,ax31,ax32]
axrow4 = [ax40,ax41,ax42]
ax = np.array([axrow0,axrow1,axrow2,axrow3,axrow4])
print(ax.shape)

# raw data
x = np.linspace(0,10,1000)
xcc = np.linspace(-x[-1]/2,x[-1]/2,len(x))
dx = x[-1]/len(x)
# offsets (delta x) between curves
x1 = 1 ; x2 = 6
yfloor = [0,0.5,1,2,4] # noise levels
for i in range(ax.shape[0]):
	print(yfloor[i])
	y1 = np.exp(-(x-x1)**2) + np.random.uniform(-yfloor[i],yfloor[i],len(x))# + yfloor[i]
	y2 = np.exp(-(x-x2)**2) + np.random.uniform(-yfloor[i],yfloor[i],len(x))# + yfloor[i]
	ax[i,0].plot(x,y1,color='r')
	ax[i,0].plot(x,y2,color='b')
	ax[i,0].locator_params(axis='y',nbins=5)
	# Cross corr
	CC = np.correlate(y1,y2,mode='same')
	ax[i,1].plot(xcc,np.flip(CC)/np.max(CC))
	ax[i,1].axvline(xcc[np.argmax(np.flip(CC))],color='k',linestyle='--')
	# SA
	SA,peak = shared_area(y1,y2,dx=dx,fitgauss=True)
	SAoff = 0
	try:
		SAneg = np.abs(SA[SA<0])
		SAoff = np.max(SAneg)
	except:
		None
	ax[i,2].plot(x,(SA+SAoff)/(np.max(SA+SAoff)))
	ax[i,2].axvline(x[np.argmax(SA)],color='k',linestyle='--')

	# formatting
	ignorey([ax[i,2]])
	if i < ax.shape[0]-1:
		ignorex([ax[i,0],ax[i,1],ax[i,2]])

ax[0,1].set_ylim(-1.2,1.2)
ax[0,1].set_xlim(-x[-1]/2,x[-1]/2)
ax[-1,0].set_xlabel(r'$t$',**tnrfont)
ax[-1,1].set_xlabel(r'$\tau$',**tnrfont)
ax[-1,2].set_xlabel(r'$\tau$',**tnrfont)
ax[2,0].set_ylabel('Data',**tnrfont)

plt.show()
fig.savefig('/home/space/phrmsf/Documents/thesis_images/shared_area/shared-area-cc_sa.png',bbox_inches='tight')
