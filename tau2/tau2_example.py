
from func_load import *



fig, axs = plt.subplots(nrows=1,ncols=3,figsize=(9,3),sharey=True)
fig.subplots_adjust(hspace=0,wspace=0)
ax = axs.ravel()

x = np.arange(0,10,0.01)
y = np.sin(x)+np.random.uniform(-1,1,len(x))
#y = np.random.random(len(x))

ax[0].scatter(x,y,color='k',marker='x',s=1)

N = 20 #xbins
M = N #200 #ybins
xbinw = (max(x)-min(x))/N
ybinw = (max(y)-min(y))/M

xbinloc = (N*x/np.max(x))//1
ybinloc = (M*y/np.max(y))//1

#for i in range(-N,N):
#	ax[1].axvline(i*xbinw,color='w',linestyle='--',linewidth=1)
#	ax[1].axhline(i*ybinw,color='w',linestyle='--',linewidth=1)

p, xedges, yedges = np.histogram2d(x,y,bins=(N,M),density=True,range=[[min(x),max(x)],[min(y),max(y)]])
im = ax[1].imshow(p.T,**kwargs,extent=[0,max(x),-2,2])

xdata = np.arange(0,10,0.1)
for off in np.arange(0,1,0.25):
	ydata = 2*np.sin(xdata+np.pi*off)
	ax[2].scatter(xdata,ydata)
	prob=0
	for i in range(len(xdata)):
		xbin = (N*(xdata[i]-np.min(x))/(np.max(x)-np.min(x)))
		ybin = (M*(ydata[i]-np.min(y))/(np.max(y)-np.min(y)))
		prob+=-2*np.log(p[int(xbin),int(ybin)])
	print(off,prob)

#ax[1].colorbar(im)
#ax[2].legend()

#xlabels = [0,2,4,6,8]
#ax[0].set_xticklabels(xlabels)
#ax[1].set_xticklabels(xlabels)
#ax[2].set_xticklabels(xlabels)
ax[0].set_ylim(-2,2)
plt.show()
