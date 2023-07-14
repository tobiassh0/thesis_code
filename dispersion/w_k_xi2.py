from list_new import *
from matplotlib.pyplot import cm

def coldpd(wc1,wp1,wc2,wp2,wce,wpe,omegas,theta=None): 
	theta = theta*const.PI/180
	sin = np.sin(theta) ; cos = np.cos(theta)
	print(theta, sin, cos)
	wpf = [wpe, wp1, wp2]
	wcf = [wce , wc1, wc2]

	l = len(omegas)
	R = np.ones(l) ; P = np.ones(l) ; L = np.ones(l) ; S = np.zeros(l) ; D = np.zeros(l) 
	B = np.zeros(l) ; F = np.zeros(l) ; A = np.zeros(l) ; C = np.zeros(l) 
	
	R = R - ((wpf[0]**2)/(omegas*(omegas + wcf[0]))) - ((wpf[1]**2)/(omegas*(omegas + wcf[1]))) - ((wpf[2]**2)/(omegas*(omegas + wcf[2])))
	L = L - ((wpf[0]**2)/(omegas*(omegas - wcf[0]))) - ((wpf[1]**2)/(omegas*(omegas - wcf[1]))) - ((wpf[2]**2)/(omegas*(omegas - wcf[2])))
	P = P -  ((wpf[0]**2)/omegas**2) -  ((wpf[1]**2)/omegas**2) - ((wpf[2]**2)/omegas**2)

	S = 0.5*(R+L) ; D = 0.5*(R-L)
	C = P*R*L
	B = R*L*(sin**2) + P*S*(1.0 +cos**2)
	F = (((R*L - P*S)**2)*(sin**4) + 4.0*(P**2)*(D**2)*(cos**2))**0.5
	A = S*(sin**2) + P*(cos**2)
	n1 = np.zeros(l, dtype=complex) ; n2 = np.zeros(l, dtype=complex) ; n3 = np.zeros(l, dtype=complex) ; n4 = np.zeros(l, dtype=complex) 
	n3 = np.lib.scimath.sqrt((R*L)/S)
	n1 =  np.lib.scimath.sqrt((B+F)/(2.0*A))
	n2 = np.lib.scimath.sqrt((B-F)/(2.0*A))
	n3 = -np.lib.scimath.sqrt((B+F)/(2.0*A))
	n4 = -np.lib.scimath.sqrt((B-F)/(2.0*A))
	del R, P, L, S, D, B, F, A
#	return n1*omegas/const.c, n2*omegas/const.c, n3*omegas/const.c, n4*omegas/const.c
	return (np.real(n1)*omegas)/const.c , (np.real(n2)*omegas)/const.c , np.real((n3*omegas)/const.c) #, (n4*omegas)/c, omegas


#i=0
#colors=['b','k']
#sim_lst = ['traceT_D_99_T_01','traceT_0_50']
#wcD = const.qe*2.1/getMass('Deuterons')
#n0 = 1e19; temp = 1000*const.qe/const.kb
#LDe = np.sqrt(const.e0*const.kb*temp/((const.qe**2)*n0))
#omegas = wcD*np.linspace(0,35,10000)
#for sim in sim_lst:
#	simloc = getSimulation('/storage/space2/phrmsf/'+sim)
#	d0 = sdfread(0)
#	#n0 = getMeanquantity(d0,'Derived_Number_Density_Electrons')
#	#nD = getMeanquantity(d0,'Derived_Number_Density_Deuterons')
#	#nT = getMeanquantity(d0,'Derived_Number_Density_Tritons')
#	_,k2,_=coldplasmadispersion(d0, 'Deuterons', 'Tritons', omegas,theta=89)
#	plt.plot(k2*LDe,omegas/wcD,color=colors[i])
##	FT2d = read_pkl('FT_2d_Magnetic_Field_Bz')
#	#d0 = sdfread(0)
#	##LDe = getLarmorRadius(d0,'Electrons')
##	klim = 0.5*2*const.PI/getdxyz(d0)
##	wlim = 0.5*2*const.PI/getdt(d0)
##	plt.imshow(np.log10(FT2d),**kwargs,cmap='magma',extent=[0,klim*LDe,0,wlim/wcD])
#	i+=1
#plt.show()


xi2 = np.linspace(0,1,8192)
#xi2 = np.array([0.01,0.5])
xi3 = 2e-3
n0 = 1e19
T0 = 1000*const.qe/const.kb
B0 = 2.1 # T
Z1 = 1 ; Z2 = 1 ; Z3 = 2
xi1 = (1/Z1)*(1-xi2*Z2-xi3*Z3)
nD = xi1*n0 ; nT = xi2*n0 ; nmin = xi3*n0
mD = getMass('Deuterons') ; mT = getMass('Tritons') ; mmin = getMass('Alphas') ; me = getMass('Electrons')
wpD = np.sqrt((nD*((Z1*const.qe)**2))/(mD*const.e0))
wpT = np.sqrt((nT*((Z2*const.qe)**2))/(mT*const.e0))
wpe = np.sqrt((n0*(const.qe**2))/(me*const.e0))
wcD = (Z1*const.qe)*B0/mD
wcT = (Z2*const.qe)*B0/mT
wce = (const.qe)*B0/me
#wcmin = (Z3*const.qe)*B0/mmin
#wpmin = np.sqrt((nmin*((Z3*const.qe)**2))/(mmin*const.e0))
LDe = np.sqrt((const.e0*const.kb*T0)/(n0*const.qe**2))
rhoM = n0*(xi1*mD+xi2*mT+xi3*mmin)
VA = B0/(const.mu0*rhoM)**.5
#plt.plot(xi2,VA/const.c) ; plt.show()
omegas = wcD*np.linspace(0,40,4096)
kMat = np.zeros((len(omegas),len(xi2)))
##plt.plot(xi2,wpD/wcD,xi2,wpT/wcD,xi2,(wpD/wpT)) ; plt.show() ; sys.exit()
for i in range(len(xi2)):
	_,k2,_=coldpd(wcD,wpD[i],wcT,wpT[i],wce,wpe,omegas,89)
	kMat[:,i] = k2

fig,ax=plt.subplots(figsize=(8,6))
#im = ax.imshow(kMat*LDe,**kwargs,cmap='jet',extent=[xi2[0],xi2[-1],omegas[0]/wcD,omegas[-1]/wcD],clim=(0,0.1))
#cbar = fig.colorbar(im)
#ax.set_ylabel(r'$\omega/\Omega_D$',fontsize=20)
#ax.set_xlabel(r'$\xi_2$',fontsize=20)
#cbar.set_label(r'$k\lambda_{De}$',fontsize=20)	
#fig.savefig('/storage/space2/phrmsf/paper/xi2_omega.png')
#fig.clf()

## plot for a fixed omega
lomega = np.arange(0,30,2)
color = cm.rainbow(np.linspace(0, 1, len(lomega)))
wcalpha = 2*const.qe*2.1/getMass('Alphas')
#plt.plot(xi2,VA/const.c,color='k')
for l, c in zip(lomega, color):
	ind = len(omegas)*l/(omegas[-1]/wcD)
	ks = kMat[int(ind),:]
	kA = (l*wcD)/VA
	kd = (ks-kA)/kA
	ksVA = ks*VA
#	ax.plot(xi2,ksVA/wcD,color=c,label=str(l//1))
#	ax.annotate(str(l//1),xy=(0.5,ksVA[len(ksVA)//2]/wcD+0.5),xycoords='data',color=c,fontsize=12)
	ax.plot(xi2,kd,color=c,label=str(l//1))
	if l>18:
		ax.annotate(str(l//1),xy=(0.5,kd[len(kd)//2]+0.025),xycoords='data',color=c,fontsize=12)

#ax.set_ylabel(r'$k_x v_A/\Omega_D$',fontsize=20)
ax.set_ylabel(r'$(k_x-k_A)/k_A$',fontsize=20)
ax.set_xlabel(r'$\xi_T$',fontsize=20)
#plt.legend(loc='best')
#fig.savefig('/storage/space2/phrmsf/paper/xi2_kperpVA.png')
fig.savefig('/storage/space2/phrmsf/paper/xi2_kdiff.png')#,bbox_inches='tight')
plt.show() ; plt.clf()




