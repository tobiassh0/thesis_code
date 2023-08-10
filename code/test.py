import numpy as np
import matplotlib.pyplot as plt 
import my_constants as const
from list_new import *
from scipy import signal
from scipy.fft import fftshift
import time
from scipy.optimize import curve_fit

def mvavg(arr,points=np.linspace(-3,3,7,dtype=int)):
	# IN 
	#	arr : 1d array
	#	points : list of points to move forward and backwards 
	#			across the data to take an average of. 
	# OUT
	#	new_arr : The moving average array
	n = len(arr)
	new_arr = np.zeros(n)
	for i in range(0,n):
		for p in points:
			if i+p > n-1 or i+p < 0:
				new_arr[i] = arr[i]
			else:
				new_arr[i] += arr[i+p]
		new_arr[i]=new_arr[i]/len(points)

	return new_arr

def tophat(length, width, offset=0, height=1):
	# IN
	#	length : the length of the top-hat array (including 0s)
	# 	width : the wdith of the top-hat for when it's >0 
	# 	offset : when to begin the top-hat in your data series (defaults to 0)
	#	height : the non-zero value of the top-hat array (defaults to 1)
	# OUT
	#	toph : the array of length with 0s and height for some width in your signal
	toph = []
	for i in range(length):
		if i/length > offset and i/length < offset+width:
			toph.append(height)
		else:
			toph.append(0)
	return toph

def mask_area(data, xlim, ylim, left, bottom, width, height, mval=0, cbar=False):
	# Mask an area within the data with a square shape starting at bottom left 
	# and moving right and upwards by width and height (in units of the axes)
	# Replaces the data array values within the shape with mval.
	# IN
	#	data : 2d array of data with shape nx, ny
	#	xlim, ylim : limits of the data in x and y (in units of the axes)
	#	left, bottom, width, height : limits of the masking area (box)
	#	mval : value to replace data within the mask
	# OUT
	#	fig : figure of the new data with an applied mask
	
	(nx, ny) = data.shape
	mask_x = nx * (left/xlim)
	mask_y = ny * (bottom/ylim)
	mask_x_width = nx * (width/xlim)
	mask_y_width = ny * (height/ylim)
	fig, axs = plt.subplot(1,2)
	ax = axs[0]
	ax.imshow(data, extent=[0,xlim,0,ylim], interpolation='nearest', origin='lower', cmap='magma')
	ax = axs[1]
	data[int(mask_x):int(mask_x+mask_x_width),int(mask_y):int(mask_y+mask_y_width)] = mval
	ax.imshow(data, extent=[0,xlim,0,ylim], interpolation='nearest', origin='lower', cmap='magma')
	if cbar:
		plt.colorbar()
	return fig

##################################
mD = const.me_to_mD
mT = const.me_to_mT

##################################

JETdata = np.loadtxt('JET26148_ICE_POWER.txt',delimiter=',')
JETpower, JETfreqs = JETdata[:,1], JETdata[:,0] # 2 columns, N rows
wnorm = 2*const.PI*17E6
print('wnorm [MHz] :: {}'.format(wnorm/(2*const.PI*1e6)))
JETfreqs = 2*const.PI*JETfreqs*1e6/(wnorm) # convert MHz to wcD
maxJETfreqs = round(max(JETfreqs))
print('MAX FREQ (JET wcD):: ',maxJETfreqs)
fig,ax=plt.subplots(figsize=(8,4))
ax.annotate(r'$\Omega_D=17$'+'MHz',xy=(0.1,0.8),xycoords='axes fraction',fontsize=20)
ax.plot(JETfreqs,JETpower,'k-')
ax.set_xlabel(r'$\omega/\Omega_D$',fontsize=20)
ax.set_ylabel('ICE intensity [dB]',fontsize=20)
ax.set_xlim(0,maxJETfreqs)
fig.savefig('JET26148_power.png',bbox_inches='tight')
fig.savefig('JET26148_power.eps',bbox_inches='tight')

#################################

#sims = ['traceT_0_00','traceT_D_99_T_01','traceT_D_89_T_11','traceT_0_50']
#print('sim','LDe','LDD','va/c','vperp/va')
#for sim in sims:
#	simloc = getSimulation('/storage/space2/phrmsf/'+sim)
#	d0 = sdfread(0)
#	Emin = 3.5e6 * const.qe
#	theta = 0.22*const.PI
#	mass = getMass('Alphas')
#	vz = np.sqrt(2*Emin/mass)
#	v_perp = vz * np.cos(theta)
#	vA = getAlfvenVel(d0)
#	LDe = getDebyeLength(d0,'Electrons')
#	LDD = getDebyeLength(d0,'Deuterons')
#	print([sim,LDe,LDD,vA/const.c,v_perp/vA])

#	wcd = getCyclotronFreq(sdfread(0),'Deuterons')
#	omegas = wcd*np.linspace(0,25,10000)
#	_,k2,_ = coldplasmadispersion(sdfread(0),'Deuterons','',omegas)
#	PIxx, PIxy, PIyy = Chi0Calc(sdfread(0),vz,k2,omegas,species='Deuterons',theta=89)
#	print(PIxx,PIxy,PIyy)

#Z1 = 1
#Z2 = Z3 = 2
#Lambda = np.arange(1e-6,1e-2,1e-7)
#Mu = np.arange(0,1,0.001)
#L,M = np.meshgrid(Lambda,Mu)
#C3 = 15
##CT = ((C3/L))*((1+1/Z1)+M*(1-Z2/Z1)+L*(1-Z3/Z1))
##plt.imshow(np.log10(CT),interpolation='nearest',origin='lower',aspect='auto',extent=[Mu[0],Mu[-1],Lambda[0],Lambda[-1]])
##plt.colorbar()
##plt.show()
#n0 = 5e19
#T0 = 2*11.6e6
#L = 3.
#dx = np.sqrt(const.kb*const.e0*T0/(n0*const.qe**2)) 
#Nx = L/dx
#CT = ((C3/Lambda))*((1+1/Z1)+0.1*(1-Z2/Z1)+Lambda*(1-Z3/Z1))
#NT = CT*Nx
#print(CT,NT)
#plt.plot(Lambda,np.log10(CT),color='k') 
#plt.plot(Lambda,np.log10(NT),color='b') 
#plt.xlabel(r'$\Lambda$')
#plt.ylabel(r'log(Particles) '+'  '+r'$[10^6]$')
#plt.show()

###############################################################################################################
###############################################################################################################




