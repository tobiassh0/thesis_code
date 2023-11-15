
#def shared_area(sig1,sig2):
#	lx = len(sig1)
#	if lx!=len(sig2): 
#		print('## ERROR ## :: Both signals need the same length')
#		raise SystemExit
#	rolled = np.zeros(lx)
#	dx = sig1[-1]/lx
#	sharea=[]
#	for i in range(lx):
#		rsig1 = np.roll(sig1,i)
#		c1 = rsig1 < sig2
#		c2 = sig2 < rsig1
#		rolled[c1] = rsig1[c1]
#		rolled[c2] = sig2[c2]
#		sharea.append(np.sum(rolled*dx))

#	return sharea

def shared_area_example():
	plt.style.use('classic')
	plt.tight_layout()
	tnrfont = {'fontsize':20,'fontfamily':'Times New Roman'}
	
	A1 = 1 ; A2 = 1
	a1 = 1 ; a2 = 6
	b1 = 5 ; b2 = 5
	x = np.linspace(0,10,1000)
	lx = len(x)
	dx = x[-1]/lx
	curve1 = A1*np.exp(-b1*(x-a1)**2)
	curve2 = A2*np.exp(-b2*(x-a2)**2)
	#plt.plot(x,curve1,x,curve2)
	
	rolled = np.zeros(lx)
	frac=50
	integral = []
	tau = [0,300,450,600]
	Areas = []
	fig = plt.figure(constrained_layout=True)
	gs = GridSpec(2,4,figure=fig)
	ax1 = fig.add_subplot(gs[0, 0])
	ax2 = fig.add_subplot(gs[0, 1])
	ax3 = fig.add_subplot(gs[0, 2])
	ax4 = fig.add_subplot(gs[0, 3])
	
	# plot stationary curve
	ax1.plot(x,curve2,color='k')
	ax2.plot(x,curve2,color='k')
	ax3.plot(x,curve2,color='k')
	ax4.plot(x,curve2,color='k')
	
	# plot moving curve
	# tau0
	rcurve1 = np.roll(curve1,tau[0])
	c1 = rcurve1 < curve2
	c2 = curve2 < rcurve1
	rolled[c1] = rcurve1[c1]
	rolled[c2] = curve2[c2]
	ax1.plot(x,rcurve1,color='b')
	ax1.fill_between(x,rolled,0,color='lightpink')
	ax1.set_xlabel(r'$t$',**tnrfont) ; ax1.set_ylabel(r'$f(t),g(t)$',**tnrfont)
	ax1.annotate(r'$\tau_1$',(0.5,0.9),xycoords='data',**tnrfont)
	Areas.append(np.sum(rolled*dx))
	#tau1
	rcurve1 = np.roll(curve1,tau[1])
	c1 = rcurve1 < curve2
	c2 = curve2 < rcurve1
	rolled[c1] = rcurve1[c1]
	rolled[c2] = curve2[c2]
	ax2.plot(x,rcurve1,color='b')
	ax2.fill_between(x,rolled,0,color='lightpink')
	ax2.set_xlabel(r'$t$',**tnrfont)
	ax2.annotate(r'$\tau_2$',(0.5,0.9),xycoords='data',**tnrfont)
	Areas.append(np.sum(rolled*dx))
	#tau2
	rcurve1 = np.roll(curve1,tau[2])
	c1 = rcurve1 < curve2
	c2 = curve2 < rcurve1
	rolled[c1] = rcurve1[c1]
	rolled[c2] = curve2[c2]
	ax3.plot(x,rcurve1,color='b')
	ax3.fill_between(x,rolled,0,color='lightpink')
	ax3.set_xlabel(r'$t$',**tnrfont)
	ax3.annotate(r'$\tau_3$',(0.5,0.9),xycoords='data',**tnrfont)
	Areas.append(np.sum(rolled*dx))
	#tau3
	rcurve1 = np.roll(curve1,tau[3])
	c1 = rcurve1 < curve2
	c2 = curve2 < rcurve1
	rolled[c1] = rcurve1[c1]
	rolled[c2] = curve2[c2]
	ax4.plot(x,rcurve1,color='b')
	ax4.fill_between(x,rolled,0,color='lightpink')
	ax4.set_xlabel(r'$t$',**tnrfont)
	ax4.annotate(r'$\tau_4$',(0.5,0.9),xycoords='data',**tnrfont)
	Areas.append(np.sum(rolled*dx))
	
	for i in range(lx):#[0,100,200,300,400,500,600,700,800,900,1000]:
		rcurve1 = np.roll(curve1,i)
		c1 = rcurve1 < curve2
		c2 = curve2 < rcurve1
		rolled[c1] = rcurve1[c1]
		rolled[c2] = curve2[c2]
		integral.append(np.sum(rolled*dx))
	#	plt.plot(x,rcurve1,x,curve2)
	#	plt.scatter(x[::frac],rolled[::frac],color='r')
	#	plt.fill_between(x,rolled,0,color='lightpink')
	#	plt.savefig('sharedarea{}.png'.format(i))
	#	plt.clf()
	ax5 = fig.add_subplot(gs[1, :])
	ax5.plot(x,integral)
	tau = np.array(tau) ; Areas = np.array(Areas)
	ax5.scatter(x[-1]*tau/lx,Areas,c='r')
	ax5.set_xlim(0,x[-1]) ; ax5.set_ylim(-0.01,1)
	ax5.set_xlabel(r'$\tau$',**tnrfont) ; ax5.set_ylabel('Area'+r'$(\tau)$',**tnrfont)
	xoff = 0.1; yoff = 0.05
	ax5.annotate(r'$A(\tau_1)$',(x[-1]*tau[0]/lx+xoff,Areas[0]+yoff),xycoords='data',color='r')
	ax5.annotate(r'$A(\tau_2)$',(x[-1]*tau[1]/lx+xoff,Areas[1]+yoff),xycoords='data',color='r')
	ax5.annotate(r'$A(\tau_3)$',(x[-1]*tau[2]/lx+xoff,Areas[2]+yoff),xycoords='data',color='r')
	ax5.annotate(r'$A(\tau_4)$',(x[-1]*tau[3]/lx+xoff,Areas[3]+yoff),xycoords='data',color='r')
	fig.savefig('shared-area.png')

	return None

if __name__=='__main__':
	from func_load import *
	from matplotlib.gridspec import GridSpec
	shared_area_example()