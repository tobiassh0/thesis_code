from func_load import *
from matplotlib.gridspec import GridSpec


sim_lst = ['lowres_D_He3/0_25_p_90']
quant = 'Magnetic_Field_Bz'

for sim in sim_lst:
	loc = getSimulation('/storage/space2/phrmsf/'+sim)
	L = getGridlen(sdfread(0))
	FT2d = read_pkl('FT_2d_'+quant)
	vA = getAlfvenVel(sdfread(0))
	times = read_pkl('times')
	dx = getdxyz(sdfread(0))
	dt = getdt(times)
	wcp = getCyclotronFreq(sdfread(0),'Protons')
	tcp = 2*const.PI/wcp
	wlim = 0.5*2*const.PI/dt
	klim = 0.5*2*const.PI/dx
	(nw,nk) = FT2d.shape
	wmax = 25*wcp ; kmax = 50*wcp/vA
	FT2d = FT2d[:int(nw*wmax/wlim),:int(nk*kmax/klim)]
	# plot
	fig, ax = plt.subplots(figsize=(8,4)) #(8,4)
	ax.imshow(np.log10(FT2d),**kwargs,extent=[0,kmax*vA/wcp,0,wmax/wcp],cmap='magma',clim=(-4,6))
	ax.set_xlabel(r'$kv_A/\Omega_p$',**tnrfont)
	ax.set_ylabel(r'$\omega/\Omega_p$',**tnrfont)
	fig.savefig('FT_2d_'+quant+'.png',bbox_inches='tight')

