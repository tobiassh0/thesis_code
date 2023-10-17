from func_load import *

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

dim = 10

X, Y = np.meshgrid([-dim, dim], [-dim, dim])
Z = np.zeros((2, 2))

angle = .5
X2, Y2 = np.meshgrid([-dim, dim], [0, dim])
Z2 = Y2 * angle
X3, Y3 = np.meshgrid([-dim, dim], [-dim, 0])
Z3 = Y3 * angle

r = 7
M = 1000
th = np.linspace(0, 2 * np.pi, M)

x, y, z = r * np.cos(th),  r * np.sin(th), angle * r * np.sin(th)

ax.plot_surface(X2, Y3, Z3, color='red', alpha=.5, linewidth=0, zorder=-1)

ax.plot(x[y < 0], y[y < 0], z[y < 0], lw=5, linestyle='--', color='red',
        zorder=0)

ax.plot_surface(X, Y, Z, color='black', alpha=.5, linewidth=0, zorder=1)

ax.plot(r * np.sin(th), r * np.cos(th), np.zeros(M), lw=5, linestyle='--',
        color='k', zorder=2)

ax.plot_surface(X2, Y2, Z2, color='red', alpha=.5, linewidth=0, zorder=3)

ax.plot(x[y > 0], y[y > 0], z[y > 0], lw=5, linestyle='--', color='red',
        zorder=4)

# plot xyz axis in centre
ax.plot([-dim-dim/10,dim+dim/10],[0,0],[0,0],color='k')
ax.plot([0,0],[-dim-dim/10,dim+dim/10],[0,0],color='k')
ax.plot([0,0],[0,0],[-dim-dim/10,dim+dim/10],color='k')
plt.axis('off')
plt.show()


#fields = ['a','b','c','d']
#for i in range(0,len(fields)):
#	ilim = len(fields)-i
#	for j in range(i,len(fields)):
#		print(i,fields[i],j,fields[j])

#simloc = getSimulation('/storage/space2/phrmsf/ECRH/ECRH_JT60U_3')
#d0 = sdfread(0)
#d1 = sdfread(1)
#
## check times
#dt = getdt(None,d0,d1)
#wce = getCyclotronFreq(d0,'Electrons')
#wcp = getCyclotronFreq(d0,'Protons')
#wpe = getPlasmaFreq(d0,'Electrons')
#wpi = getPlasmaFreq(d0,'Protons')
#tce = 2*const.PI/wce
#tcp = 2*const.PI/wcp
#times = read_pkl('times')
#print(dt/tce,dt/tcp)
#
## check freq lim
#wUH = getUpperHybrid(wce,wpe)
#wlim = 0.5*2*const.PI/dt
#dw = 2*wlim/len(times)
#print(wlim/wce,wlim/wcp,wlim/wUH)
#print(dw/wce,dw/wcp)
