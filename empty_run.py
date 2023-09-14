from func_load import *
from matplotlib.gridspec import GridSpec


sim_lst = ['lowres_D_He3/0_25_p_90']
quant = 'Magnetic_Field_Bz'

for sim in sim_lst:
	loc = getSimulation('/home/space/phrmsf/Documents/EPOCH/epoch-4.17.16/epoch1d/old/0006qp')#/storage/space2/phrmsf/'+sim)
	d0 = sdfread(0)
	wcD = const.qe*2.1/getMass('Deuterons')
	omegas = wcD*np.linspace(0,25,1000)
	k1, k2, k3 = coldplasmadispersion(d0, 'Deuterons', 'Tritons', omegas, theta=89)


	L = getGridlen(d0)
	FT2d = read_pkl('FT_2d_'+quant)
	vA = getAlfvenVel(d0)
	times = read_pkl('times')
	dx = getdxyz(d0)
	dt = getdt(times)
	wcp = getCyclotronFreq(d0,'Deuterons')
	tcp = 2*const.PI/wcp
	wlim = 0.5*2*const.PI/dt
	klim = 0.5*2*const.PI/dx
	(nw,nk) = FT2d.shape
	wmax = wlim#25*wcp
	kmax = klim#50*wcp/vA
	dw = 0.5*2*const.PI/times[-1]
	dk = 0.5*2*const.PI/L
	print(dw/getCyclotronFreq(d0,'Deuterons'))
#	FT2d = FT2d[:int(nw*wmax/wlim),:int(nk*kmax/klim)]
	# plot
	fig, ax = plt.subplots(figsize=(8,4)) #(8,4)
	ax.imshow(np.log10(FT2d),**kwargs,extent=[0,kmax*vA/wcp,0,wmax/wcp],cmap='magma',clim=(-4,6))
	ax.set_xlabel(getWavenumberLabel('Deuterons'),**tnrfont)
	ax.set_ylabel(r'$\omega/$'+getOmegaLabel('Deuterons'),**tnrfont)
	plt.show()
	#fig.savefig('FT_2d_'+quant+'.png',bbox_inches='tight')

