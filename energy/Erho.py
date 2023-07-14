import sdf
import numpy as np
from list_new import *#
from list_new import *
import matplotlib.pyplot as plt
import my_constants as const

sim_loc = getSimulation('/storage/space2/phrmsf/traceT_0_01')
ind_lst = list_sdf(sim_loc)
Ex = getQuantity1d(sdfread(0),'Electric_Field_Ex')
nt = len(ind_lst) ; nx = Ex.shape[0]

### technique 1
E_rho = np.zeros((nt,nx))
rho_c = getMeanquantity(sdfread(0),'Derived_Number_Density_Electrons') * const.qe
dx = getdxyz(sdfread(0))
r = dx/2

def ELECTRICEX(index):
	print(index)
	return getQuantity1d(sdfread(index),'Electric_Field_Ex')

import multiprocessing as mp

pool = mp.Pool(mp.cpu_count()-2)
#Bz = np.array(pool.map_async(ELECTRICEX,ind_lst).get(99999))
Ex = np.array(pool.map_async(ELECTRICEX,ind_lst).get(99999))
#dEx = Ex - np.roll(Ex,1)
#E_rho = (rho_c/const.e0 - dEx/dx)*(r/2)
pool.close()
plt.imshow((Ex),interpolation='nearest',origin='lower',aspect='auto',cmap='seismic',extent=[0,nx,0,nt])
plt.colorbar() ; plt.show()

#B0 = getMeanField3D(sdfread(0),'Magnetic_Field_B')
#Ey = np.zeros((nt,nx)); Ez = np.zeros((nt,nx))
#By = np.zeros((nt,nx))
#Bz = np.zeros((nt,nx))
#Ex = np.zeros((nt,nx))
#sval=[]
#for j in range(0,nt):
#	if j%1000 == 0: print(j)
#	dj = sdfread(j)
#	rho_c = getQuantity1d(dj,'Derived_Charge_Density')
#	ex = getQuantity1d(dj,'Electric_Field_Ex')
#	E_rho[j,:] = ((rho_c/const.e0) - (ex/dx))*dx/2
#	By[j,:] = getQuantity1d(dj,'Magnetic_Field_By')
#	Bz[j,:] = getQuantity1d(dj,'Magnetic_Field_Bz')
#	Ex[j,:] = getQuantity1d(dj,'Electric_Field_Ex')
#	try:
#		Ey[j,:] = getQuantity1d(dj,'Electric_Field_Ey')
#		Ez[j,:] = getQuantity1d(dj,'Electric_Field_Ez')
#		sval.append(j)
#	except:
#		None

##plt.imshow(np.log10(E_rho),interpolation='nearest',origin='lower',aspect='auto',cmap='jet',extent=[0,nx,0,nt])
##plt.colorbar()
##plt.savefig('Erho.png')

##plt.clf()
#Ecollect= np.sqrt(Ey**2 + Ez**2)
#Ecollect = Ecollect[~np.all(Ecollect==0,axis=1)] # remove lines with 0 in them (including the first row?)
##plt.imshow(np.log10(Ecollect),interpolation='nearest',origin='lower',aspect='auto',cmap='jet',extent=[0,nx,0,nt])
##plt.colorbar()
##plt.savefig('Ecollect.png')

#Ecr = []
#diff = sval[1] - sval[0]
#dumpfiles(sval,'sval')
#sval = np.arange(0,len(Ecollect),1)
#for i in range(len(sval)):
#	Ecr.append(E_rho[diff*i,:]/Ecollect[i,:]) 

#dumpfiles(E_rho,'E_rho')

##plt.clf()
##plt.imshow(np.log10(Ecr),interpolation='nearest',origin='lower',aspect='auto',cmap='jet',extent=[0,nx,0,nt])
##plt.colorbar()
##plt.savefig('Ecollect_Erho.png')
##print(np.mean(Ecr),np.std(Ecr))

#Ey = Ey[~np.all(Ey==0,axis=1)] # remove zero lines
#Ez = Ez[~np.all(Ez==0,axis=1)]
#r = np.mean(abs(Ey)/abs(Ez))
#a = np.sqrt(1/(1+1/(r**2)))
#b = a/r
#Bx = np.sqrt(B0**2 - Bz**2 - By**2)
#Eyr = a*E_rho
#Ezr = b*E_rho
#print('ratios',r,a,b)
#Sx = (1/const.mu0)*(Eyr*Bz - Ezr*By)
#Sy = -1*(1/const.mu0)*(Ex*Bz-Ezr*Bx)
#Sz = (1/const.mu0)*(Ex*By - Eyr*Bx)
#Sx = Sx[1:,:]
#Sy = Sy[1:,:]
#Sz = Sz[1:,:]
#SmagInfer = np.sqrt(Sx*Sx + Sy*Sy + Sz*Sz) # magnitude of the Poynting Vector
#plt.imshow(np.log10(SmagInfer),aspect='auto',interpolation='nearest',origin='lower')
#plt.show()
#plt.savefig('SmagInfer.png')



#for j in range(len(sval)):	
#	print(str(np.around(100*j/len(sval),2))+'%')
#	d = sdfread(int(diff*sval[j]))
#	Ex = getQuantity1d(d,'Electric_Field_Ex')
#	Ey = getQuantity1d(d,'Electric_Field_Ey')
#	Ez = getQuantity1d(d,'Electric_Field_Ez')
#	Bx = getQuantity1d(d,'Magnetic_Field_Bx')
#	By = getQuantity1d(d,'Magnetic_Field_By')
#	Bz = getQuantity1d(d,'Magnetic_Field_Bz')
#	Sx[j,:] = (1/const.mu0)*(Ey*Bz - Ez*By)
#	Sy[j,:] = -1*(1/const.mu0)*(Ex*Bz-Ez*Bx)
#	Sz[j,:] = (1/const.mu0)*(Ex*By - Ey*Bx)
## remove first line
#Sx = Sx[1:,:]
#Sy = Sy[1:,:]
#Sz = Sz[1:,:]
#Smag = np.sqrt(Sx*Sx + Sy*Sy + Sz*Sz) # magnitude of the Poynting Vector

#FTSmag = Poynting2dFT(times,nt,nx,plot=False)
#plt.imshow(np.log10(FTSmag),aspect='auto',interpolation='nearest',origin='lower')
#plt.show()


## technique 2
#calc = True
#while calc:
#	for i in range(1,nt):
#		try:
#			Ex0 = getQuantity1d(sdfread(i),'Electric_Field_Ex')
#			Jx0 = getQuantity1d(sdfread(i),'Current_Jx')
#			Ey0 = getQuantity1d(sdfread(i),'Electric_Field_Ey')
#			Jy0 = getQuantity1d(sdfread(i),'Current_Jy')
#			Ez0 = getQuantity1d(sdfread(i),'Electric_Field_Ez')
#			Jz0 = getQuantity1d(sdfread(i),'Current_Jz')
#			print('file found :: ',i)
#			calc = False
#			break
#		except:
#			calc = True
#sigma_x = (1/nx)*(Ex0/Jx0) 
#sigma_y = (1/nx)*(Ey0/Jy0) 
#sigma_z = (1/nx)*(Ez0/Jz0) 
#plt.plot(np.linspace(0,nx,nx),sigma_x,color='b')
#plt.plot(np.linspace(0,nx,nx),sigma_y,color='g')
#plt.plot(np.linspace(0,nx,nx),sigma_z,color='r')
#plt.show()
#print(np.mean(sigma_x),np.std(sigma_x))
#print(np.mean(sigma_y),np.std(sigma_y))
#print(np.mean(sigma_z),np.std(sigma_z))
#sigma_x = np.mean(sigma_x)
#sigma_y = np.mean(sigma_y)
#sigma_z = np.mean(sigma_z)
#sigma = np.sqrt(sigma_x**2 + sigma_y**2 + sigma_z**2) # pythagorean summation/contribution

#species = getIonSpecies(sdfread(0))

#Ex = np.zeros((nt,nx))
#for j in range(0,nt):
#	if j%10 == 0: print(j)
#	dj = sdfread(j)
#	rho_c = getQuantity1d(dj,'Derived_Charge_Density')
#	v=0
#	for spec in species:
#		try:
#			EnAvg = getQuantity1d(dj,'Derived_Average_Particle_Energy_'+spec)
#			v += np.sqrt(2*EnAvg/getMass(spec)) # TODO; fix this as I dont think it's just a summation of velocities (might be?)
#		except:
#			None
#	Ex[j,:] = (rho_c*v)/(sigma)

#plt.imshow(np.log10(Ex),interpolation='nearest',origin='lower',aspect='auto',cmap='jet')
#plt.colorbar()
#plt.show()


## parallelized
#sim_loc = getSimulation('/storage/space2/phrmsf/cold_JET26148')
#times = read_pkl('times')
#nt = times.shape[0]
#print(nt)
#Ex0 = getQuantity1d(sdfread(0),'Electric_Field_Ex')
#nx = Ex0.shape[0]
#print(nx)

#calc = True
#while calc:
#	for i in range(1,nt):
#		try:
#			Ex0 = getQuantity1d(sdfread(i),'Electric_Field_Ex')
#			Jx0 = getQuantity1d(sdfread(i),'Current_Jx')
#			Ey0 = getQuantity1d(sdfread(i),'Electric_Field_Ey')
#			Jy0 = getQuantity1d(sdfread(i),'Current_Jy')
#			Ez0 = getQuantity1d(sdfread(i),'Electric_Field_Ez')
#			Jz0 = getQuantity1d(sdfread(i),'Current_Jz')
#			print('file found :: ',i)
#			calc = False
#			break
#		except:
#			calc = True
#sigma_x = (1/nx)*(Ex0/Jx0) 
#sigma_y = (1/nx)*(Ey0/Jy0) 
#sigma_z = (1/nx)*(Ez0/Jz0) 
##	print(np.mean(sigma_x),np.std(sigma_x))
##	print(np.mean(sigma_y),np.std(sigma_y))
##	print(np.mean(sigma_z),np.std(sigma_z))
#sigma_x = np.mean(sigma_x)
#sigma_y = np.mean(sigma_y)
#sigma_z = np.mean(sigma_z)
#sigma = np.sqrt(sigma_x**2 + sigma_y**2 + sigma_z**2) # pythagorean summation/contribution
#	
#def func(i):
#	if i%10 == 0: print(i)
#	di = sdfread(i)
#	rho_c = getQuantity1d(di,'Derived_Charge_Density')
#	v = 0
#	species = ['Deuterons','Tritons']
#	for spec in species:
#		EnAvg = getQuantity1d(di,'Derived_Average_Particle_Energy_'+spec)
#		v += np.sqrt(2*EnAvg/getMass(spec)) # TODO; fix this as I dont think it's just a summation of velocities (might be?)
#	return (rho_c*v)/(sigma)

#nt_arr = np.arange(0,nt,1)
#num_cores = multiprocessing.cpu_count()
#print('num cores :: ',num_cores)
##Ex = np.zeros((nt,nx))
#Ex = Parallel(n_jobs=num_cores)(delayed(func)(i) for i in nt_arr)
#print(Ex)
#plt.imshow(np.log10(Ex),interpolation='nearest',origin='lower',aspect='auto',cmap='jet')
#plt.show()




















