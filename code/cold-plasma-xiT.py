from func_load import *
import time

def now():
	return time.strftime("%Y%m%d-%H%M%S")

def getnewest():
	oslst = os.listdir()
	lst1 = [i for i in oslst if 'nsol1' in i]	
	lst2 = [i for i in oslst if 'nsol2' in i]	
	return lst1[0].replace('.pkl',''), lst2[0].replace('.pkl','')
	

def colddisp(B0=2.1,n0=5e19,xialpha=1e-3,wmax=25,dtheta=89):
	zd = 1
	zt = 1
	ze = -1
	zalpha = 2
	md = const.me_to_mD
	mt = const.me_to_mT
	malpha = const.me_to_malpha
	wcD = zd*const.qe*B0/(const.me*md)
	wcT = zt*const.qe*B0/(const.me*mt)
	wcalpha = zalpha*const.qe*B0/(malpha*const.me)
	wce = abs(ze*const.qe*B0/(const.me))
	wpe = np.sqrt((n0*const.qe**2)/(const.me*const.e0))
	wpalpha = np.sqrt((n0*xialpha*(const.qe*zalpha)**2)/(const.me*malpha*const.e0))

	omegas = wcD*np.linspace(0,wmax,100)
	theta = dtheta*const.PI/180
	xiT = np.linspace(0,1,100)
	FREQ, XIT = np.meshgrid(omegas,xiT)

	zeta = (zd/md)*(1-zalpha*xialpha)
	sin2 = (np.sin(theta))**2
	cos2 = (np.cos(theta))**2

	## R terms
	cr1 = -(wpe**2/FREQ)*((zt**2)/(mt*(FREQ+wcT))-(zt*zd**2)/(md*(FREQ+wcD)))
	cr0 = 1-(1/FREQ)*((wpe**2)/(FREQ+wce)+(wpalpha**2)/(FREQ+wcalpha)+zeta*(wpe**2)/(FREQ+wcD))
	R = cr1*XIT + cr0

	## L terms
	cl1 = -(wpe**2/FREQ)*((zt**2)/(mt*(FREQ-wcT))-(zt*zd**2)/(md*(FREQ-wcD)))
	cl0 = 1-(1/FREQ)*((wpe**2)/(FREQ-wce)+(wpalpha**2)/(FREQ-wcalpha)+zeta*(wpe**2)/(FREQ-wcD))
	L = cl1*XIT + cl0

	## P terms
	cp1 = -(wpe**2/FREQ**2)*((zt**2/mt)-((zt*zd**2)/md))
	cp0 = 1-(1/FREQ**2)*(wpalpha**2+(wpe**2)*(1+zeta))
	P = cp1*XIT + cp0

	## S terms
	cs1 = -(zt*wpe**2/(2*FREQ))*((zt/mt)*((2*FREQ)/(FREQ**2-wcT**2))-(zd**2/md)*((2*FREQ)/(FREQ**2-wcD**2)))
	cs0 = 1-(1/(2*FREQ))*((wpalpha**2)*(2*FREQ/(FREQ**2-wcalpha**2))+2*(FREQ*wpe**2)*((1/(FREQ**2-wpe**2))+zeta*(1/(FREQ**2-wcD**2))))
	S = cs1*XIT + cs0

	## D terms
	cd1 = ((zt*wpe**2)/(2*FREQ))*((zt/mt)*(2*wcT/(FREQ**2-wcT**2))-(zd**2/md)*(2*wcD/(FREQ**2-wcD**2)))
	cd0 = (1/(2*FREQ))*((wpalpha**2)*(2*wcalpha/(FREQ**2-wcalpha**2))+2*(wpe**2)*(wce/(FREQ**2-wce**2)+zeta*(wcD/(FREQ**2-wcD**2))))
	D = cd1*XIT + cd0

	## A terms
	ca1 = cs1*sin2 + cp1*cos2
	ca0 = cs0*sin2 + cp0*cos2
	A = ca1*XIT + ca0

	## B terms
	cb2 = cr1*cl1*sin2 + cp1*cs1*(1+cos2)
	cb1 = (cr1*cl0+cr0*cl1)*sin2 + (cp1*cs0+cp0*cs1)*(1+cos2)
	cb0 = cr0*cl0*sin2 + cp0*cs0*(1+cos2)
	B = cb2*XIT**2 + cb1*XIT + cb0

	## F terms
	cf4 = ((cr1*cl1-cp1*cs1)**2)*(sin2**2) + cos2*(2*cp1*cd1)**2
	cf3_1 = 2*(cr1*cl1-cp1*cs1)*((cr1*cl0-cr0*cl1)-(cp1*cs0+cp0*cs1))*(sin2**2)
	cf3_2 = 8*(cd1*cd0*(cp1**2)+cp1*cp0*(cd1**2))*cos2
	cf3 = cf3_1+cf3_2
	cf2_1 = (2*(cr0*cl0-cp0*cs0)*(cr1*cl1-cp1*cs1)+((cr1*cl0-cr0*cl1)-(cp1*cs0+cp0*cs1))**2)*(sin2**2)
	cf2_2 = 4*((cp1**2)*(cd0**2)+(cd1**2)*(cp0**2)+4*cp1*cp0*cd1*cd0)*cos2
	cf2 = cf2_1+cf2_2
	cf1_1 = (2*((cr1*cl0-cr0*cl1)-(cp1*cs0+cp0*cs1))*(cr0*cl0-cp0*cs0))*(sin2**2)
	cf1_2 = 8*(cp1*cp0*(cd0**2)+cd1*cd0*(cp0**2))*cos2
	cf1 = cf1_1+cf1_2
	cf0 = ((cr0*cl0-cp0*cs0)**2)*(sin2**2)+cos2*(2*cp0*cd0)**2
	F2 = (cf4*XIT**4)+(cf3*XIT**3)+(cf2*XIT**2)+cf1*XIT+cf0
	F = np.sqrt(F2)

	nsol1 = np.real(np.lib.scimath.sqrt((B+F)/(2*A)))
	nsol2 = np.real(np.lib.scimath.sqrt((B-F)/(2*A)))
	tnow = now()
	dumpfiles(nsol1,'nsol1-'+tnow)
	dumpfiles(nsol2,'nsol2-'+tnow)
	return nsol1, nsol2

if __name__ == '__main__':
	try:
		nfile1, nfile2 = getnewest()
		nsol1 = read_pkl(nfile1)
		nsol2 = read_pkl(nfile2)
	except:
		nsol1, nsol2 = colddisp()
	plt.imshow(nsol2,**kwargs)
	plt.show()











