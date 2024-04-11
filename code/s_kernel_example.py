from func_load import * 
import matplotlib.gridspec as gridspec

#setup meshgrid
lim = 10
x = np.linspace(0,lim,500)
y = x
X,Y = np.meshgrid(x,y)
extents = [0,lim,0,lim]

# 90 deg
#Z0 = Y
Z90 = 2*X
#plt.imshow(Z90,**kwargs,cmap='Greys',extent=extents) ; plt.show()

# 45 deg
Z45 = 2*(X + Y)
#plt.imshow(Z45,**kwargs,cmap='Greys') ; plt.show()

# sinusoidal at 45 deg
Zsin = np.sin((X+Y))
##plt.imshow(Zsin,**kwargs,cmap='Greys',extent=extents)
#kmag,kang = Kernel(Zsin,'sobel')
#kang = kang*180/const.PI # deg
#plt.imshow(kang[1:-1,1:-1],**kwargs,cmap='Greys',clim=(-360,360))
#cbar = plt.colorbar()
#plt.show() ; sys.exit()

Zarr = [Z90,Z45,Zsin]

# plot figure
fig = plt.figure()
width = 12
rows = len(Zarr)
columns = 3 # (Z, sobel, scharr)
# setup gridspec
gs1 = gridspec.GridSpec(rows, width//columns, left=0.05, right=0.40, wspace=0., hspace=0.)# right spacing creates space for gs2 
gs2 = gridspec.GridSpec(rows, 2*width//columns, left=0.53, right=0.98, wspace=0., hspace=0.) # nrows, ncols, l, r, wsp, hsp

# Z plots
zlabels = [r'$90^\circ$',r'$135^\circ$',r'$\sin(x+y)$']
axZ = []
for i in range(0,rows):
	axi = fig.add_subplot(gs1[i,int(0.3*width//columns):]) # [row number, :3]
	axZ.append(axi)

for z in range(len(axZ)):
	# plot Z heatmaps
	axZ[z].imshow(Zarr[z],**kwargs,extent=extents,cmap='Greys')
	# formatting ~~~
	axZ[z].text(1,7,zlabels[z],va="bottom",ha="left",fontsize=16,
				bbox=dict(facecolor='white',edgecolor='black',boxstyle='square,pad=0.15'))
	axZ[z].set_ylabel(r'$y$',**tnrfont)
	if z == len(axZ)-1:
		axZ[z].set_yticks([0,2,4,6,8,10])
		axZ[z].set_xticklabels([0,2,4,6,8,10])
		axZ[z].set_xlabel(r'$x$',**tnrfont)
	else:
		axZ[z].set_yticks([2,4,6,8,10])
		axZ[z].set_xticklabels([])
	# ~~~

# kang hist plots
axsobel = [] ; axscharr = []
for i in range(0,rows):
	axso = fig.add_subplot(gs2[i,:width//columns]) # [row number, 3:6]
	axsc = fig.add_subplot(gs2[i,width//columns:])   # [row number, 6:]
	axsobel.append(axso)
	axscharr.append(axsc)

for k in range(len(axsobel)):
	# calc kernel gradients
	_,ksobel = Kernel(Zarr[k],'sobel')
	ksobel = ksobel*180/const.PI # deg
	axsobel[k].hist(ksobel.flatten(),bins=360,range=(-180,180),density=True)
	_,kscharr = Kernel(Zarr[k],'sobel')
	kscharr = kscharr*180/const.PI # deg
	# plot hist of kernel grads (angles)
	axscharr[k].hist(kscharr.flatten(),bins=360,range=(-180,180),density=True)
	# formatting ~~~
	axsobel[k].set_xlim(-180,180)
	axscharr[k].set_xlim(-180,180)
	axsobel[k].set_ylim(0,1.25)
	axscharr[k].set_ylim(0,1.25)
	#
	axsobel[k].set_xticks([-180,-90,0,90,180])
	axsobel[k].set_xticklabels([])
	axscharr[k].set_xticks([-180,-90,0,90,180])
	axscharr[k].set_xticklabels([])
	# 
	axsobel[k].set_yticks([0,0.25,0.5,0.75,1])
	axsobel[k].set_yticklabels([0,0.25,0.5,0.75,1])
	axscharr[k].set_yticks([0,0.25,0.5,0.75,1])
	axscharr[k].set_yticklabels([])
	if k == 2:
		axsobel[k].set_xticks([-180,-90,0,90,180])
		axsobel[k].set_xticklabels([-180,-90,0,90,180])
		axsobel[k].set_xlabel(r'$\theta$'+'  '+'[deg]',**tnrfont)
		axsobel[k].set_yticks([0,0.25,0.5,0.75,1])
		axsobel[k].set_yticklabels([0,0.25,0.5,0.75,1])
		#
		axscharr[k].set_xticks([-90,0,90,180])
		axscharr[k].set_xticklabels([-90,0,90,180])
		axscharr[k].set_xlabel(r'$\theta$'+'  '+'[deg]',**tnrfont)
	# ~~~
axsobel[0].set_title('Sobel',**tnrfont)
axsobel[1].set_ylabel('Normalised count',**tnrfont)
axscharr[0].set_title('Scharr',**tnrfont)
# plt.show()
fig.savefig('/home/space/phrmsf/Documents/thesis_images/example_kernels.png',bbox_inches='tight')
