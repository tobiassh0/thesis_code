
from func_load import *


fig, axs = plt.subplots(nrows=2,ncols=2,figsize=(10,8),layout='constrained')
# fig.subplots_adjust(hspace=0,wspace=0)
ax = axs.ravel()

x = np.arange(0,10,0.005)
y = np.sin(x)
for i in range(len(x)):
    y[i] += np.random.normal(y[i],1)
#y = np.sin(x)+np.random.normal(-1,1,len(x))
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
print('SUM HERE ', np.sum(p))
# normalise as A.A.Breimann 2021
floor = 1/((max(y)-min(y))*(max(x)-min(x))) # area of x and y plane
weight = 0.7
p = ((1-weight)*floor + weight*p)/np.sum(p)
print('SUM HERE ', np.sum(p))
im = ax[1].imshow(p.T,**kwargs,extent=[0,max(x),-2,2])

xdata = np.arange(0,10,0.1)
phase = np.arange(-1,1.25,0.25)
tau2 = np.zeros(len(phase))
for j in range(len(phase)):
	ydata = 2*np.sin(xdata+np.pi*phase[j])
	if phase[j] < 0:
		color = 'r'
	elif phase[j] > 0:
		color = 'b'
	else:
		color = 'k'
	ax[2].plot(xdata,ydata,color=color,alpha=(1-abs(phase[j]))**1.5)
	prob=0
	for i in range(len(xdata)):
		xbin = (N*(xdata[i]-np.min(x))/(np.max(x)-np.min(x)))-1
		ybin = (M*(ydata[i]-np.min(y))/(np.max(y)-np.min(y)))-1
		prob+=-2*np.log(p[int(xbin),int(ybin)])
	print(phase[j],prob)
	tau2[j] = prob

ax[3].plot(phase,tau2,'-o')

#xlabels = [0,2,4,6,8]
#ax[0].set_xticklabels(xlabels)
#ax[1].set_xticklabels(xlabels)
#ax[2].set_xticklabels(xlabels)

# # shared y-axes 
# ax[0].set_ylim(-2,2)

# x and y labels
ax[0].set_ylabel(r'$y$',**tnrfont)
ax[0].set_xlabel(r'$x$',**tnrfont)
ax[1].set_xlabel(r'$x$',**tnrfont)
ax[1].set_ylabel(r'$y$',**tnrfont)
ax[2].set_xlabel(r'$x$',**tnrfont)
ax[2].set_ylabel(r'$y$',**tnrfont)
ax[3].set_xlabel(r'$x_{off}$',**tnrfont)
ax[3].set_ylabel(r'$\tau^2$',**tnrfont)

# panel labels
label = ['(a)','(b)','(c)','(d)']
for j in range(len(ax)):
	ax[j].annotate(label[j],xy=(0.05,0.9),xycoords='axes fraction',**tnrfont)

if phase[np.argmin(tau2)] == 0:
	fig.savefig('tau2_sinusoid_example.png')

plt.show()
