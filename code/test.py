
from func_load import *
import correlation.cross_correlation as cc

x = np.linspace(0,10,100)
sig1 = np.exp(-(x-1)**2)
sig2 = np.exp(-(x-5)**2)
croc = cc.getCrossCorrelation(sig1,sig2,name='test')
fig,ax = cc.plotCrossCorrelation(np.linspace(-x[-1]/2,x[-1]/2,len(x)),croc)
plt.show()

