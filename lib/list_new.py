''' 
	Still need to import the functions here too otherwise list_new doesn't work.
	Could try and fix this, could create a package out of this whole directory. Haven't
	done it yet, may get around to it.
'''

import os, sys
import numpy as np
import sdf
import pickle
import matplotlib.pyplot as plt
from scipy import special as spec
from scipy import stats, signal
from scipy.optimize import curve_fit
from matplotlib.path import Path
import multiprocessing as mp ## parallelisation
#from cv2 import cv2 # used to make images

## OLD PACKAGES USED
#import scipy.fftpack
#from itertools import cycle
#from scipy.optimize import curve_fit
#from matplotlib.lines import Line2D 
#from scipy.ndimage.filters import gaussian_filter

## 
from bispectral_analysis import *
import my_constants as const
## 

#plt.rcParams['text.usetex'] = True
#plt.rcParams['text.latex.unicode'] = True
plt.style.use('classic')
plt.tight_layout()
plt.rcParams['axes.formatter.useoffset'] = False

def formatting():
	kwargs ={'interpolation':'nearest','origin':'lower','aspect':'auto'}
	tnrfont = {'fontsize':20,'fontname':'Times New Roman'}
	return kwargs, tnrfont

kwargs, tnrfont = formatting()
#---------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------#

# Returns the file name/loc of the simulation being analysed
def getSimulation(loc=''):
	cwd = os.getcwd()
	if loc == '':
		sim_file = input('Input Simulation Directory:\n!>>')
	else:
		sim_file = loc
	try:
		os.chdir(sim_file)
	except:
		print('# ERROR # : Simulation file does not exist.')
		raise SystemExit
	return os.path.join(cwd,sim_file)

## returns all of the ion species in a simulation from looking at file0, could also do a function where user inputs each 
# of these but would be annoying to type each time
def getIonSpecies(file0):
	keys = getKeys(file0) # file0.__dict__.keys() # returns all the names of quantities in the simulation sdf file
	# check something which is always dumped like derived charge dens
	names = []
	for k in keys:
		if 'Derived_Number_Density_' in k and 'Electrons' not in k:
			names.append(k)
		else:
			continue
	names = [s.replace('Derived_Number_Density_','') for s in names]
	len_names = len(names)
	if len_names <= 2:
		maj2spec = ''
		minspec = ''
		if len_names == 1:
			return names[0], maj2spec, minspec # case when there is only one ion species
		elif len_names == 2:
			return names[0], maj2spec, names[1] # case when there are two ions total (maj & min)
		else:
			print('# ERROR # : Can\'t count number of species in sim.')
			raise SystemExit
	elif len_names == 3:
		return names[0], names[1], names[2] # maj1, maj2, min # i.e. D, T, alphas
	else:
		print('# ERROR # : More than three ion species in sim.')
		raise SystemExit

def getAllSpecies(file0):
	keys = getKeys(file0) # file0.__dict__.keys() # returns all the names of quantities in the simulation sdf file
	# check something which is always dumped like derived charge dens
	names = []
	for k in keys:
		if 'Derived_Number_Density_' in k:
			names.append(k)
		else:
			continue
	names = [s.replace('Derived_Number_Density_','') for s in names]
	return names #all species in sim, in order of how they appear in EPOCH output (usually electrons, maj1, maj2, min)

# Scans and returns a list of files as objs, readable in the form "files[0], files[1] , ..." 
def filelist(lim):
	# creates array of sdf read objects (files) either for 0000.sdf or 00000.sdf file labelling
	try:	
		files = [sdf.read(('%04d'%i)+'.sdf') for i in range(0,lim)]
	except:
		try:
			files = [sdf.read(('%05d'%i)+'.sdf') for i in range(0,lim)]
		except:
			files = np.asarray(None) # None array
	return np.array(files) 

# Lists all the sdf files in the given simulation, converts this to a list of indices and returns it (as list of numbers)
def list_sdf(sim_file_loc):
	sdf_list = [i for i in os.listdir(sim_file_loc) if ".sdf" in i and 'TempOut' not in i] # ignores Temp files 
	index_list = np.sort(np.asarray([s.replace('.sdf','') for s in sdf_list], dtype=int))
	return index_list

# Reads a given sdf file as per the index
def sdfread(index,d=0,l5=True,l4=False):
	try:
		d=sdf.read(('%05d'%index)+'.sdf')
		l5=True ; l4=not l5
	except:
		try:
			d=sdf.read(('%04d'%index)+'.sdf')	
			l4=True ; l5=not l4
		except:
			if index==0:
				index+=1 # load 1st file instead of 0th 
				try:
					if l4:
						d=sdf.read(('%04d'%index)+'.sdf')			
					else:
						d=sdf.read(('%05d'%(index))+'.sdf')
				except:
					print('\# ERROR \# : Can\'t load {}.sdf'.format(index))
					d=None
	return d

def getKeys(d):
	keys = d.__dict__.keys() # returns all the names of quantities in the simulation sdf file
	return keys

# Gives info about the simulation grid (used for dxdydz and more later)
def getGrid(d):
	return d.__dict__["Grid_Grid_mid"].data

# Takes in the grid shape and size then returns the total length across x,y and z (depending on dimensionality)
def getGridlen(d):
	# TODO: Make more general so that we can define lx,ly,lz and then return them all based on dimensions?
	l = getGrid(d)
	if len(l) == 1:
		lx = l[0]
		del l
		return lx[-1] - lx[0]
	elif len(l)==2:
		lx = l[0] ; ly = l[1]
		del l
		return lx[-1]-lx[0] , ly[-1]-ly[0]
	elif len(l)==3:
		lx = l[0] ; ly = l[1] ; lz = l[2]
		del l
		return lx[-1]-lx[0] , ly[-1]-ly[0], lz[-1]-lz[0]

# Returns the physical size of domain steps in x,y,z (depending on dimension of sim)
def getdxyz(d):
	l=getGrid(d)
	if (len(l)==1):
		lx = l[0]
		dx = (lx[-1]-lx[0])/(len(lx)-1)
		del lx ,l
		return dx
	elif (len(l)==2):
		lx = l[0]
		dx = (lx[-1]-lx[0])/(len(lx)-1)
		ly = l[1]
		dy = (ly[-1]-ly[0])/(len(ly)-1)
		del lx, ly , l
		return dx, dy
	elif (len(l)==3):
		lx = l[0]
		dx = (lx[-1]-lx[0])/(len(lx)-1)
		ly = l[1]
		dy = (ly[-1]-ly[0])/(len(ly)-1)
		lz = l[2]
		dz = (lz[-1]-lz[0])/(len(lz)-1)
		del lx, ly, lz , l
		return dx, dy, dz

# Returns a list of times (seconds) from each sdf file
def getTimes(files):
	# Scans all sdf objs in the "files", returns the times and stores into an array
	times=np.zeros((len(files)))
	for i in range(0,len(files)):
		times[i] = files[i].__dict__['Header']['time']
		# print(files[i].__dict__['Header']['time'])
	# times[0] =float(0) #if you start from t=0, first timestep be small but not zero, so make it
			#0. Remove if not starting from 0. There is probably a better way of doing this
	# could re-code this to not use the file list
	return times

# returns the mean dt value across a simulation, assuming equal time-steps 
def getdt(times=None,file0=None,file1=None):
	try:
		dt = (times[-1]-times[0])/len(times)
	except:
		try:
			dt = (file1.__dict__['Header']['time']) - (file0.__dict__['Header']['time'])
		except:
			dt = None
	return dt

# plots an array of dt values through time
def plotdt(times=None):
	if not bool(times.any()):
		times = read_pkl('times')
	dt = times - np.roll(times,1)
	plt.plot(times[:-2],dt[:-2]) ; plt.ylabel('dt',**tnrfont) ; plt.xlabel('time',**tnrfont) 
	plt.axhline(np.mean(dt[1:-10]),linestyle='--',color='k')	
	print(np.mean(dt[1:-10]))
	plt.show()
	plt.savefig('dt_times.png')
	plt.clf()
	return None

# appends an array of times batch-wise
def batch_getTimes(times,start,stop):
	if start == stop:
		file = sdfread(len(index_list)-1) 
		times[-1] = file.__dict__['Header']['time']
	for i in range(start,stop+1):
		file = sdfread(i)
		times[i] = file.__dict__['Header']['time']
		# print(files[i].__dict__['Header']['time'])
	# times[0] =float(0) #if you start from t=0, first timestep be small but not zero, so make it
			#0. Remove if not starting from 0. There is probably a better way of doing this
	# could re-code this to not use the file list, also depends on index_list which is not passed but presumably a variable
	del file
	return times

# Gets the total simulation ranges for a len(filelist) < 1024
def getSimRanges(files):
	grid_length = getGridlen(files[0])
	N_grid_points = len(getGrid(files[-1])[0])
	times = getTimes(files)
	duration = times[-1] - times[0]
	N_time_points = len(files)
	return grid_length, N_grid_points, duration, N_time_points

# Gets the simulation ranges from a batch-wise loaded array
def batch_getSimRanges(SIM_DATA):
	index_list,file0,filelast,times = SIM_DATA
	grid_length = getGridlen(file0)
	N_grid_points = len(getGrid(filelast)[0])
	duration = times[-1] - times[0]
	N_time_points = len(index_list)
	return grid_length, N_grid_points, duration, N_time_points

# Gets dispersion limits of a simulation for a len(filelist) < 1024
def getDispersionlimits(files):
	grid_length, N_grid_points, duration, N_time_points = getSimRanges(files)
	wfrac = 2.0*const.PI/duration 		#2pi/time of sim		(2pi/T)	
	kfrac = 2.0*const.PI/grid_length 		#2pi/length of grid 	(2pi/L)
	klim = kfrac*0.5*N_grid_points
	wlim = wfrac*0.5*N_time_points
	# highest freq we can resolve # see OneNote lab-book for more info (utilises Nyquist theorem)
	return  klim, wlim #,kfrac, wfrac

# Gets dispersion limits of a simulation for a batch-wise loaded array
def batch_getDispersionlimits(SIM_DATA):
	grid_length, N_grid_points, duration, N_time_points = batch_getSimRanges(SIM_DATA)
	wfrac = 2.0*const.PI/duration 		# delta omega # 2pi/T	
	kfrac = 2.0*const.PI/grid_length 		#2pi/length of grid 	(2pi/L)
	klim = kfrac*0.5*N_grid_points
	wlim = wfrac*0.5*N_time_points
	# highest freq we can resolve # see OneNote lab-book for more info (utilises Nyquist theorem)
	del grid_length, N_grid_points, duration, N_time_points, wfrac, kfrac
	return  klim, wlim

def norm_DispersionLimits(DISP_DATA):
	index_list, file0, filelast, times, klim, wlim, wc_maj, wce, va, lambdaD, wpe, wpi, wnorm = DISP_DATA
	grid_length, N_grid_points, duration, N_time_points = batch_getSimRanges((index_list,file0,filelast,times))
	# klim,wlim = batch_getDispersionlimits((self.index_list,self.file0,self.filelast,self.times)) # non-normalised units
	_,_,klim_prime,wlim_prime,tlim_prime = norm_PlottingLimits((index_list,file0,filelast,times),klim,wlim,wc_maj,wce,va,lambdaD,wpe,wpi,wnorm)
	## k,w,t_prime are now normalised according to whether there is a alfven wave present (hence B field)
	## TODO; change this so that I can choose the normalisation
	return klim_prime, wlim_prime, tlim_prime

# Returns the quantity of a given time-step
def getQuantity1d(d, quantity):
	quan = d.__dict__[quantity]
	# quan_name, quan_units = quan.name, quan.units
	# x = quan.grid.data[0]
	# x = np.array(x)
	quan_array = np.array(quan.data)
	return quan_array

# Returns the mean of a quantity
def getMeanquantity(d,quantity): 
	mean_quan = getQuantity1d(d,quantity)
	mean_quan = np.mean(mean_quan)
	return mean_quan

# Returns the mean of a given field: e.g. quantity = 'Electric_Field_E' or 'Magnetic_Field_B'
def getMeanField3D(d,quantity): #include all of the header apart from last charchter (x,y,z) which is the coord
	try:
		mean_x = getMeanquantity(d,quantity+'x')
	except:
		mean_x = 0
	
	try:
		mean_y = getMeanquantity(d,quantity+'y')
	except:
		mean_y = 0

	try:
		mean_z = getMeanquantity(d,quantity+'z')
	except:
		mean_z = 0

	#doesn't account for sign so will return lower even if amplitude is large and negative
	# modified so that it will try to read the component, returning 0 otherwise
	return (mean_x**2+mean_y**2+mean_z**2)**0.5

# Gets the mass (m) of a species defined by a dictionary
def getMass(species):
	masses = {'Electrons': const.me,
	'Electron': const.me,
	'FElectrons': const.me,
	'Left_Electrons': const.me,
	'Right_Electrons': const.me, 
	'Protons': const.me*const.me_to_mp, 
	'Proton': const.me*const.me_to_mp, 
	'PProtons': const.me*const.me_to_mp, 
	'Alphas': const.me*const.me_to_malpha, 
	'Alpha': const.me*const.me_to_malpha, 
	'Deuterons': const.me*const.me_to_mD,
	'Deuteron': const.me*const.me_to_mD, # sometimes singular
	'Deutrons': const.me*const.me_to_mD, # sometimes misspelled
	'Tritium' : const.me*const.me_to_mT,
	'Tritons' : const.me*const.me_to_mT,
	'Triton' : const.me*const.me_to_mT,
	'Helium3': const.me*const.me_to_He3, 
	'He3': const.me*const.me_to_He3,
	'Borons': const.me*const.me_to_B11,
	'Boron': const.me*const.me_to_B11,
	'B11': const.me*const.me_to_B11,
	'Ions': const.me*const.me_to_mp} # assuming by Ions the user means Protons

	if species not in masses:
		if species == '':
			return 0
		else:
			print('Species [{}] mass is not in dictionary, check name passed for spelling mistakes'.format(species))
		raise SystemExit
	else: 
		return masses.get(species)

# Gets the charge number (Z) of a species defined by a dictionary
def getChargeNum(species):
	charges = {'Electrons': 1,
	'FElectrons': 1,
	'Left_Electrons': 1,
	'Right_Electrons': 1, 
	'Protons': 1, 
	'PProtons': 1, 
	'Alphas': 2, 
	'Alpha': 2,
	'Deuterons': 1,
	'Deutrons': 1, # sometimes misspelled
	'Tritium': 1,
	'Tritons': 1,
	'Helium3': 2, 
	'He3': 2,
	'Borons': 5,
	'Boron': 5,
	'B11': 5,
	'Ions': 1,
	'': 0}

	if species not in charges:
		if species == '':
			print('No species provided, Z=0')
		else:
			print('Species [{}] charge number is not in dictionary, check name passed for spelling mistakes'.format(species))
		raise SystemExit
	else: 
		return charges.get(species)

def getIonlabel(species):
	labels = {'Electrons': r'$e$',
	'FElectrons': r'$e$', # fast electrons
	'Left_Electrons': r'$e$',
	'Right_Electrons': r'$e$', 
	'Protons': r'$p$', 
	'PProtons': r'$p$', # sometimes have two of the same species
	'Alphas': r'$\alpha$', 
	'Alpha': r'$\alpha$', 	
	'Helium4': r'$\alpha$',
	'He4': r'$\alpha$',
	'Deuterons': r'$D_2$',
	'Deutrons': r'$D_2$', # sometimes misspelled
	'Tritium': r'$T_3$',
	'Tritons': r'$T_3$',
	'Helium3': r'$He_3$', 
	'He3': r'$He_3$',
	'Borons': r'$B_{11}$',
	'Boron': r'$B_{11}$',
	'B11': r'$B_{11}$',
	'Ions': 'Ions'}
	
	if species not in labels:
		if species == '':
			print('No maj2 species, no Ion label provided')
		else:
			print('Species [{}] label is not in dictionary, check name passed for spelling mistakes'.format(species))
			raise SystemExit
	else: 
		return labels.get(species)
	
def getOmegaLabel(species):
	labels = {'Electrons': r'$\Omega_e$',
	'FElectrons': r'$\Omega_e$', # fast electrons
	'Left_Electrons': r'$\Omega_e$',
	'Right_Electrons': r'$\Omega_e$', 
	'Protons': r'$\Omega_p$', 
	'PProtons': r'$\Omega_p$', 
	'Alphas': r'$\Omega_\alpha$', 
	'Alpha': r'$\Omega_\alpha$', 	
	'Deuterons': r'$\Omega_D$',
	'Deutrons': r'$\Omega_D$', # sometimes misspelled
	'Tritium': r'$\Omega_T$',
	'Tritons': r'$\Omega_T$',
	'Helium3': r'$\Omega_{He3}$', 
	'He3': r'$\Omega_{He}$',
	'Borons': r'$\Omega_{B11}$',
	'Boron': r'$\Omega_{B11}$',
	'B11': r'$\Omega_{B11}$',
	'Ions': r'$\Omega_i$',
	'Ion': r'$\Omega_i$'}
	
	if species not in labels:
		if species == '':
			print('No maj2 species, no Omega label provided')
		else:
			print('Species [{}] label is not in dictionary, check name passed for spelling mistakes'.format(species))
			raise SystemExit
	else: 
		return labels.get(species)

def getWavenumberLabel(species):
	return r'$v_A/$'+getOmegaLabel(species)

	
def getEffectiveMass(d):
	ions = list(filter(None,getIonSpecies(d))) # remove empty elements
	narr = np.zeros(len(ions)) ; marr = np.zeros(len(ions))
	for i in range(len(ions)):
		narr[i] = getMeanquantity(d,'Derived_Number_Density_'+ions[i])
		marr[i] = getMass(ions[i])

	n0 = getMeanquantity(d,'Derived_Number_Density_Electrons')
	meff = np.sum(marr*narr/n0) ## (xi1*m1 + xi2*m2 + xi3*m3)
	return meff

def getEffectiveCyclotronFreq(d):
	meff = getEffectiveMass(d)
	return const.qe*getMeanField3D(d,'Magnetic_Field_B')/meff
	
def getEffectivePlasmaFreq(d):
	meff = getEffectiveMass(d)
	n0 = getMeanquantity(d,'Derived_Number_Density_Electrons')
	return np.sqrt((n0*(const.qe)**2)/(meff*const.e0))

# Gets the cyclotron frequency of a species. Is generalised so can work for positive and negatively charged species
def getCyclotronFreq(d,species,Z=None):
	if not Z: Z = getChargeNum(species)
	return Z*const.qe*getMeanField3D(d,"Magnetic_Field_B")/getMass(species)

# Returns the plasma freq of a given species (generalised to take Derived Number Density, will fail without this in first dump)
def getPlasmaFreq(d,species):
	Z = getChargeNum(species)
	return ((getMeanquantity(d,'Derived_Number_Density_'+species)*((Z*const.qe)**2))/(getMass(species)*const.e0))**0.5

# Calculates thermal velocity, species and mass are hardcoded to an electron.
def getThermalVelocity(d,species):
	temp = getTemperature(species)
	return (const.kb*temp/getMass(species))**0.5

# Returns the Alfven Velocity of a simulation assuming one species (given) dominates
def getAlfvenVel(d):
	spec_lst = getIonSpecies(d)
	totden = 0
	for spec in spec_lst:
		try:
			totden += getMass(spec)*getMeanquantity(d,'Derived_Number_Density_'+spec) # loops through all ionic species in the sim including minority
		except:
			continue
	return getMeanField3D(d,'Magnetic_Field_B')/(const.mu0*totden)**0.5

# Returns the Debye Length at a given time-step (careful to use this in loops as temperature and density will change)
def getDebyeLength(d, species): # Species should always in general be Electrons
	temp = getTemperature(species)
	return (const.e0*const.kb*temp/(getMeanquantity(d,"Derived_Number_Density_"+species)*const.qe**2))**0.5 

def getTemperature(species):
	temp = 0
	try:
		TempSDF = sdf.read('TempOut00000.sdf')
		temp = getMeanquantity(TempSDF,'Derived_Temperature_'+species) # should read as K
		return temp #/(const.eV_to_K*1E3) ## returns temp in Kelvin
	except:
		try:
			temp = getMeanquantity(sdfread(0),'Derived_Temperature_'+species) 
			# could loop through files but temperature changes
#			print('Temperature of {} [keV] :: {}'.format(species,str(temp/(const.eV_to_K*1E3))))
			return (temp)
		except:
			if temp == 0:
#				print('# ERROR # : Can\'t get temperature from files...')
				temp = 1
#				temp = float(input('Temperature of '+species+'? [keV] :: '))
				return temp*const.eV_to_K*1E3  #in keV
			else:
				print('# ERROR # :: TEMP')# shouldnt have to get here
	
# Returns the Larmor radius of a species depending on the mean magnetic field, and their average energy
def getLarmorRadius(d, species):
#	v_perp, v_para = getPerpParaVel(d,species) # TODO; finish this
#	return getMass(species)*v_perp/(getChargeNum(species) * getMeanField3D(d,'Magnetic_Field_B')) # should be more accurate than just using energy and assuming even split
	try:
		ek = getMeanquantity(d, 'Derived_EkBar_'+species)		
	except:
		ek = getMeanquantity(d, 'Derived_Average_Particle_Energy_'+species)
	return np.sqrt(2*getMass(species)*ek)/(const.qe*getChargeNum(species)*getMeanField3D(d,'Magnetic_Field_B'))
	
# Rotate velocities so they align with the Magnetic Field
def RotateVel(d, species, phi=None):
	if not phi: phi_x, phi_y = getMagneticAngle(d)
	mass = getMass(species)
	vx = getQuantity1d(d,'Particles_Px_'+species)/mass
	vy = getQuantity1d(d,'Particles_Py_'+species)/mass
	vz = getQuantity1d(d,'Particles_Pz_'+species)/mass

	vx_prime = vx*np.cos(phi_x) - vz*np.sin(phi_x)
	vy_prime = vy
	vz_prime = vx*np.sin(phi_x) + vz*np.cos(phi_x)
	
	return vx_prime, vy_prime, vz_prime

# Get the perp and para components to velocity of a given species
def getPerpParaVel(d, species):
	try:
		bx = getMeanField3D(d,'Magnetic_Field_Bx')
		by = getMeanField3D(d,'Magnetic_Field_By')
		bz = getMeanField3D(d,'Magnetic_Field_Bz')
	except:
		print('# ERROR # : Magnetic fields not readable') ; return None, None
	mass = getMass(species)
	try:
		vx = getMeanquantity(d,'Particles_Px_'+species)/mass
		vy = getMeanquantity(d,'Particles_Py_'+species)/mass
		vz = getMeanquantity(d,'Particles_Pz_'+species)/mass
		E0 = getMeanquantity(d,'Derived_Average_Particle_Energy_'+species)
		vbool = v2 == (vx**2 + vy**2 + vz**2)
		print(vbool)
		v2 = 2*E0/mass
	except:
		print('# ERROR # : Particle momentum/energy not readable') ; return None, None
	vpara = (vx*bx+vy*by+vz*bz)/(bx**2 + by**2 + bz**2)
	vperp = np.sqrt(v2 - vpara)
	return vpara, vperp

#	minE = getQuantity1d()
#	phi_x, phi_y = getMagneticAngle(d)
#	try:
#		vy = getMeanquantity(d, 'Particles_Py_'+species)/getMass(species)
#		vx = getMeanquantity(d, 'Particles_Px_'+species)/getMass(species)
#	except:
#		print('Using energy for Larmor...') # have to assume energy is split evenly between x and y (wont be)
#		try:
#			ek = getMeanquantity(d, 'Derived_EkBar_'+species)		
#		except:
#			ek = getMeanquantity(d, 'Derived_Average_Particle_Energy_'+species)
#		vx = vy = (1/np.sqrt(2)) * np.sqrt(2*ek/getMass(species)) # TODO ; finish this calculation including vz components
#	vperp_x = vx * np.sin(phi_x); vpara_x = np.cos(phi_x)
#	vperp_y = vy * np.sin(phi_y); vpara_y = np.cos(phi_y)
#	
#	v_perp = np.sqrt(vperp_x**2 + vperp_y**2) ; v_para = np.sqrt(vpara_x**2 + vpara_y**2)
#	return v_perp, v_para

# Returns the plotting limits of a 2d FT fieldmatrix which looks to see if there is a magnetic field present (va flag) or not
def PlottingLimits(files, species, Z):
	wci   = getCyclotronFreq(files[0],species,Z)
	va   = getAlfvenVel(files[0])
	lambD= getDebyeLength(files[0],'Electrons')
	wpe	 = getPlasmaFreq(files[0],species='Electrons')
	wpp  = getPlasmaFreq(files[0],species=species)
	# wp   = np.sqrt(wpp*2+wpe*2)
	# u0x  = getThermalVelocity(files[0],'Electrons')# abs(getMeanquantity(files[0],'Particles_Vx_Electrons')) # TODO: Change this form being hardcoded to just that of Electrons
	# print('Check here:',u0x, wpe, u0x/wpe, lambD)
	
	if wci or va != 0:
		print('Using alfven & wci normalisation...')
		klim = getDispersionlimits(files)[0]*(va/wci)
		wlim = getDispersionlimits(files)[1]*(1/wci)
		tlim = getSimRanges(files)[2]*(wci/(2*const.PI))
	else:
		print('Using lambda & wpe normalisation...')
		klim = getDispersionlimits(files)[0]*lambD
		wlim = getDispersionlimits(files)[1]*(1/wpe)
		tlim = getSimRanges(files)[2]*(wpe/(2*const.PI))
	# klim = getDispersionlimits(files)[0]*(u0x/wpe)
	# klim = getDispersionlimits(files)[0]*(va/wci)
	# wlim = getDispersionlimits(files)[1]*(1/wci)
	# tlim = getSimRanges(files)[2]*(wci/(2*const.PI))

	print('duration: ', getSimRanges(files)[2])
	print('va=',va,'wci=',wci,'wpe=',wpe)
	print('klim=',klim,'wlim=',wlim,'tlim=',tlim)
	return wci, va, klim, wlim, tlim

# Normalises plotting limits to the klim and wlim for a batch-wise loaded array
def norm_PlottingLimits(SIM_DATA,klim,wlim,wci,wce,va,lambD,wpe,wpp,wnorm):
	#index_list, file0, filelast, times = SIM_DATA
	if wci or va != 0:
		print('Using alfven & wci normalisation...')
		klim_prime = batch_getDispersionlimits(SIM_DATA)[0]*(va/wci)
		wlim_prime = batch_getDispersionlimits(SIM_DATA)[1]*(1/wnorm)
		tlim_prime = batch_getSimRanges(SIM_DATA)[2]*(wnorm/(2*const.PI))
	else:
		print('Using lambda & wpe normalisation...')
		klim_prime = batch_getDispersionlimits(SIM_DATA)[0]*lambD
		wlim_prime = batch_getDispersionlimits(SIM_DATA)[1]*(1/wnorm)
		tlim_prime = batch_getSimRanges(SIM_DATA)[2]*(wnorm/(2*const.PI))
	# klim = getDispersionlimits(files)[0]*(u0x/wpe)
	# klim = getDispersionlimits(files)[0]*(va/wci)
	# wlim = getDispersionlimits(files)[1]*(1/wci)
	# tlim = getSimRanges(files)[2]*(wci/(2*const.PI))

	print('va=',va,'wci=',wci,'wpe=',wpe)
	print('klim_prime=',klim_prime,'wlim_prime=',wlim_prime,'tlim_prime=',tlim_prime)
	return wci, va, klim_prime, wlim_prime, tlim_prime


# Returns a matrix of t vs x of a quantity. Can be used to get an EM field before an FT and/or visualise an entire sim in one plot
def getfieldmatrix(files,quantity): 
	fieldmatrix = np.zeros([len(files),len(getGrid(files[0])[0])]) # note it relies on another function
	for i in range(0,len(files)): #loop over all files/timesteps
		fieldmatrix[i][:] = getQuantity1d(files[i],quantity)	#populates each row with the field values at that time
	return fieldmatrix #returns matrix with field values at each t,x


# Plots a variable "varname" through x-space at a single time-step 
def plot1d(d,varname):
	fig,ax=plt.subplots(figsize=(5,5))
	var=d.__dict__[varname] #pull data we want
	x=var.grid.data[0] #without the '[0]' this is a tuple with only 1 element
	ax.plot(x[1:], var.data, 'r+-') #staggered grid, the length has one more value than all quantities found so far, throw 1st away...
	ax.set_xlabel(var.grid.labels[0] + r' $(' + var.grid.units[0] + ')$',fontsize='12') #labels are in the sdf file, may need to reconsider x units
	ax.set_ylabel(var.name + r' $(' + var.units + ')$',fontsize='12') #as are the units]
	# fig.savefig(var.name+'.jpeg')
	fig.savefig(home_path+varname+'.jpeg',bbox_inches='tight')
	ax.clear()

# Moving average calculation of a quantity
def temporal_average(index_list,quantity,points=np.linspace(-1,1,3,dtype=int)):
	# TODO : Make this use the file_list rather than index_list and then phase this index_list out from the whole code
	# Moving average calculation of a quantity over N-points. Evenly weighted in time (backwards and forwards) and discards points outside of range of
	# index list. e.g. points=[-1,0,1] will carry out moving average either side of data point (including itself '0').
	valtspec=[]
	time=[]
	n = max(index_list)
	for index in index_list:
		# time.append(getQuantity1d(sdfread(index),'Wall_time'))
		vald = 0.
		counter = 0.
		for p in points:
			if index+p < 0 or index+p > n:
				continue
			else:
				dp = sdfread(index+p)
				vald += abs(getMeanquantity(dp,quantity))
				counter += 1
				del dp
		valtspec.append(vald/counter)
	del n, vald

	return valtspec, time


# Plots the velocity vs. real-space diagram (Phase space) using the list of files found in setup(). Then saves them to home dir  
def plotPhaseSpace(files,species=False):
	print('PLOTTING PHASE SPACE')
	d0 = files[0]
	if not species:
		species = getAllSpecies(d0)
	print(len(species))
	color = ['r','b','g']
	times = getTimes(files)
	LambdaD = getDebyeLength(d0,species[0])
	print('LambdaD = ',LambdaD)
	w_p = getPlasmaFreq(d0,species[0])
	Omega_cp = getCyclotronFreq(d0,species[0])
	times = (times*w_p)/(2*const.PI)
	v_th = getThermalVelocity(d0,species[0])
	del d0
	# setup figure
	fig, ax = plt.subplots(figsize=(10,5))
	# loop through files
	xmin = 0 ; xmax = 0
	for i in range(len(files)):
		print(i)
		for j in range(len(species)):
			vx_name = files[i].__dict__['Particles_Vx_'+species[j]]
			x_name = vx_name.grid.data[0]/LambdaD 		# Normalises distance to units of Debye length
			ax.scatter(x_name,vx_name.data/v_th,color=color[j],s=1.) 		# Normalises speed to v_thermal at t=0
			if np.min(x_name) < xmin:
				xmin = np.min(x_name)
			if np.max(x_name) > xmax:
				xmax = np.max(x_name)
		ax.set_xlabel(r'$x/\lambda_D$',fontsize=18)
		ax.set_ylabel(r'$v_x/v_{th}$', fontsize=18)
		ax.set_yticks([-180,-135,-90,-45,0,45,90,135,180])
		ax.set_ylim(-180,180)
		ax.set_ylim(-120,120)
		ax.set_xlim(xmin,xmax)
#		ax.set_xticks([0,500,1000,1500,2000,2500])
#		ax.set_xticklabels(['0','500','1000','1500','2000','2500'])
		ax.set_title(r'$t\omega_{pe}/2\pi = $'+str(np.around(times[i],3)), fontsize=18)
		fig.savefig('velocity{}.jpeg'.format(i), bbox_inches="tight")
		ax.clear()

# Plots the electric field potential in the x direction using the numpy.gradient() function so as to maintain array size.
def plotPhi(files,species='Electrons'):
	print('PLOTTING PHI THROUGH TIME')	
	d0 = files[0]
	LambdaD = getDebyeLength(d0,species=species)
	w_p = getPlasmaFreq(d0)
	Omega_cp = getCyclotronFreq(d0,'Protons',1)
	times = getTimes(files)*(w_p)/(2*const.PI)
	del d0

	fig,ax=plt.subplots(figsize=(10,5))
	i=0
	for file in files:
		EX = file.__dict__['Electric_Field_Ex']
		x = EX.grid.data[0]
		dx = getdxyz(file)
		# phi = -1*np.gradient(EX, dx)
		EX = np.array(EX.data)
		phi = -1*EX*dx
		
		ax.plot(x[1:]/LambdaD,phi)

		ax.set_ylabel(r'$\phi$' +' ' +r'$ [V]$',fontsize=18)
		ax.set_xlabel(r'$x/\lambda_D$',fontsize=18)
		# ax.set_ylim(-200,200)
		# ax.set_yticks([-7.5,-6.0,-4.5,-3.0,-1.5,0,1.5,3.0,4.5,6.0,7.5])
		# ax.set_yticklabels([r'$-7.5$',r'$-6.0$',r'$-4.5$',r'$-3.0$',r'$-1.5$',r'$0.0$',r'$1.5$',r'$3.0$',r'$4.5$',r'$6.0$',r'$7.5$'])
		ax.set_title(r'$t\omega_{pe}/2\pi = $'+str(np.around(times[i],3)), fontsize=18)
		
		fig.savefig(home_path+'phi{}.jpeg'.format(i), bbox_inches='tight')
		ax.clear()
		print(i)
		i+=1

# Plots the mean of the absolute electric field at each time step. From this we can plot the instability growth rate (characterised by gamma)  
def plotAveEX(files,species):
	print('PLOTTING MEAN OF ABSOLUTE EX')
	Etot=[]
	d0 = files[0]
	w_p = getPlasmaFreq(d0)
	gamma = w_p/2
	LambdaD = getDebyeLength(d0,species)
	print('gamma: '+str(gamma))
	OFFSET = 15.7
	times = getTimes(files)
	dt = times[1] - times[0]
	GAMMA = np.exp((gamma*times)-OFFSET)
	# times = (times*w_p)/(2*const.PI)
	del d0
	
	for i in range(len(files)): 
		EX = getQuantity1d(files[i],'Electric_Field_Ex')
		EXmean = np.mean(np.abs(EX))
		Etot.append(EXmean)

	EXunits = files[0].__dict__['Electric_Field_Ex'].units
	fig, ax = plt.subplots(figsize=(6,6))

	## Finding optimum offset
	# first = np.where(times>0.01)[0][0]
	# last = np.where(times<0.1)[-1][-1]
	# print(first,last)
	# Etot_dummy = Etot[first:last]
	# time_dummy = np.linspace(0.025,0.1,len(Etot_dummy))
	# pars, cov = curve_fit(f=exponential, xdata=time_dummy, ydata=Etot_dummy, p0=[0], bounds=(-np.inf, np.inf))
	# print(pars)
	# ax.plot(time_dummy,exponential(time_dummy,pars[0]),linestyle='--',color='b')

	ax.plot(times,Etot,color='b',linestyle='-',alpha=1.)
	# ax.plot(times[first[0][0]:last[-1][-1]],GAMMA[first[0][0]:last[-1][-1]],linestyle='--',color='r')
	ax.plot(times,GAMMA,linestyle='--',color='r')

	left, bottom, width, height = [0.41, 0.22, 0.45, 0.2]
	ax2 = fig.add_axes([left, bottom, width, height])
	ax2.plot(times,Etot,color='b',linestyle='-')
	ax2.plot(times,GAMMA,linestyle='--',color='r')
	ax2.set_xlim(0,0.08)
	ax2.set_ylim(1E-6,1E-3)
	ax2.set_yscale('log')
	ax2.tick_params(axis='both', which='major', labelsize=10)
	# ax2.set_ylabel(r'$\langle|E_x|\rangle $'+'  '+r'$ [ $'+ EXunits +r'$ ]$',fontsize=14)
	# ax2.set_xlabel(r'$t$' + ' ' +r'$ [s] $',fontsize=14)


	ax.set_ylim(1E-6,1E-2)
	ax.set_yscale('log')
	ax.set_ylabel(r'$\langle|E_x|\rangle $'+'  '+r'$ [ $'+ EXunits +r'$ ]$',fontsize=18)
	ax.set_xlabel(r'$t$' + ' ' +r'$ [s] $',fontsize=18)
	ax.annotate(r'$\gamma = {}$'.format(np.around(gamma,1)), xy=(0.25,0.85), xycoords='axes fraction', color='r')

	fig.savefig(home_path+'EXmean_t.jpeg',bbox_inches='tight')

# Plots any general quantity that covers time and space then plots a normalised heat map of this value (if norm='on') and saves it accordingly.
def plotNormHeatmap(matrix,xlabel='x',ylabel='y',cbar_label='',extent=[None,None],xylim=((None,None),(None,None)),norm=1,cbar=True,cmap='jet'):
	# normalise
	matrix = matrix/norm
	
	# check limits
	if extent.count(None) != 0:
		extent = matrix.shape # matrix shape is extent of plot
	
	if xylim.count(None) != 0:
		xylim=((0,matrix.shape[0]),(0,matrix.shape[1]))

	# plot
	fig,ax=plt.subplots(figsize=(8,6))
	im = plt.imshow(matrix,**kwargs,extent=extent,cmap=cmap)

	# limits and labels
	if cbar:
		tcbar = plt.colorbar(label=cbar_label)
		for y in tcbar.ax.get_yticklabels():
			y.set_fontsize(12)
	ax.set_xlim(xylim[0][0],xylim[0][1])
	ax.set_ylim(xylim[1][0],xylim[1][1])
	ax.set_xlabel(xlabel,**tnrfont)
	ax.set_ylabel(ylabel,**tnrfont)

	return fig,ax



# =================================================
# =================================================

def HanningWindowT(fieldmatrix):
	# look at the "np.hanning documentation"
	new_field = np.zeros((fieldmatrix.shape[0],fieldmatrix.shape[1]))
	han = np.hanning(fieldmatrix.shape[0]) #window with length of x dim
	for i in range(0, fieldmatrix.shape[1]): #loop over t vignal.spectrogramalues
		new_field[:,i] = fieldmatrix[:,i]*han
	del han
	return new_field

def HanningWindow2D(field):
	#new_field = np.zeros((field.shape[0],field.shape[1]))
	hant = np.hanning(field.shape[0]) #window with length of x dimHanningWindowT
	hanx = np.hanning(field.shape[1])
	han2D = np.outer(hant,hanx)
	new_field = field*han2D
	return new_field

def HanningWindowK(field):
	new_field = np.zeros((field.shape[0],field.shape[1]))
	han = np.hanning(field.shape[1]) #window with length of x dim
	for i in range(0, field.shape[0]): #loop over t values
		new_field[i,:] = field[i,:]*han
	return new_field

def BartWindowT(fieldmatrix):
	new_field = np.zeros((fieldmatrix.shape[0],fieldmatrix.shape[1]))
	han = np.bartlett(fieldmatrix.shape[0]) #window with length of x dim
	for i in range(0, fieldmatrix.shape[1]): #loop over t vignal.spectrogramalues
		new_field[:,i] = fieldmatrix[:,i]*han
	del han
	return new_field

def BlackWindowT(fieldmatrix):
	new_field = np.zeros((fieldmatrix.shape[0],fieldmatrix.shape[1]))
	han = np.blackman(fieldmatrix.shape[0]) #window with length of x dim
	for i in range(0, fieldmatrix.shape[1]): #loop over t vignal.spectrogramalues
		new_field[:,i] = fieldmatrix[:,i]*han
	del han
	return new_field
	
def HammWindowT(fieldmatrix):
	new_field = np.zeros((fieldmatrix.shape[0],fieldmatrix.shape[1]))
	han = np.hamming(fieldmatrix.shape[0]) #window with length of x dim
	for i in range(0, fieldmatrix.shape[1]): #loop over t vignal.spectrogramalues
		new_field[:,i] = fieldmatrix[:,i]*han
	del han
	return new_field


def Windows(field,window):
	win_dict = {
	"bartlett" : BartWindowT,
	"blackman" : BlackWindowT,
	"hamming" : HammWindowT,
	"hanning" : HanningWindowT
	}
	return win_dict[window](field)  # win_dict.get(window,lambda : "# ERROR # : Invalid window name.\n")
	
## TODO; make a new function where the user can choose what type of normalisation they want (by default this is just read in atm using va as a trigger)
# def Norm_Type():

# Get the 1d FT of field data
def get1dTransform(fieldmatrix, window=False, start=0): 
	# takes log of absolute shifted 1D FT (x -> k) of the matrix from getfieldmatrix()
	if window: fieldmatrix = HanningWindowT(fieldmatrix)
	else: print('!# Warning #! : 1d FFT is not being Hanning windowed.')

	preshiftFT = np.zeros((fieldmatrix.shape[0], fieldmatrix.shape[1]),complex)
	for t in range(0, fieldmatrix.shape[0]):
		preshiftFT[t][:]=np.fft.fft(fieldmatrix[t][:])
	
	shift = np.fft.fftshift(preshiftFT[start:,start:],1)
	FT = np.abs(shift) # Destroy phase information
	FT_half = FT[:,int(FT.shape[1]/2):] # All time, k>0
	del preshiftFT, FT
	return FT_half


# Plot the 1d FT of field data (using the above get1d function) and plots it with time on the y-axis
def plot1dTransform(FT_matrix,klim,tlim,klabel=r'$v_A/\Omega_D$',wlabel=r'$\Omega_i$',cbar=False,cmap='seismic',clim=(-3,4.5)): # Plots t vs k as a heat map of a field quantity
	# Pass it the matrix from "get1dTransform". Also pass it klim and tlim from "PlottingLimits()"
	trFT = np.log10(FT_matrix[:][1:])
	extent=[0, klim, 0, tlim]
	fig, ax = plt.subplots(figsize=(8,8))
	im = plt.imshow(trFT,**kwargs,extent=extent,cmap=cmap,clim=clim)
	
	ax.set_xlabel(r'$k$'+klabel,**tnrfont) #need to change when normalising 
	ax.set_ylabel(r'$t$'+wlabel+r'$/2\pi$',**tnrfont)

	if (cbar): plt.colorbar()
	del trFT

	return fig, ax


# Get the 2d FT of field data
def get2dTransform(fieldmatrix,window=True):
	if window: preshift = np.fft.fft2(HanningWindowT(fieldmatrix))[:,:] # windowed in T by default
	else: preshift = np.fft.fft2(fieldmatrix)[:,:] 
	shift = np.fft.fftshift(preshift) # shift to zero-freq so is symmetrical (read numpy fft.fftshift documentation)
	shiftchopped = shift[int(fieldmatrix.shape[0]/2):,int(fieldmatrix.shape[1]/2):] # only take positive frequencies (w) and wavenumbers (k)
	# shiftchopped = shift[:,:] # plots the whole FFT space
	shiftchopped = np.abs(shiftchopped) # destroy phase information
	del preshift, shift
	return shiftchopped

# Plot the 2d FT of field data using the above get2d function
def plot2dTransform(FFT_matrix,klim,wlim,klabel=r'$v_A/\Omega_D$',wlabel=r'$\Omega_i$',cbar=False,clim=(-4,6),cmap='magma'):
	# In:
	#	FFT_matrix , 2d FFT matrix of a field quantity e.g. Magnetic_Field_Bz
	#	klim , the extent to want to plot and show the heatmap # normalised
	#	wlim , same as above but for the frequency component # normalised
	#	cbar , colour bar boolean
	# Out:
	#	returns the figure and axis plotted so can add own annotations etc
	tr = np.log10(FFT_matrix)[1:,1:]
	extent=[0,klim,0,wlim]
	fig, ax = plt.subplots(figsize=(8,4)) #(8,4)
	ax.set_xlabel(r'$k$'+klabel,**tnrfont) #need to change when normalising 
	ax.set_ylabel(r'$\omega/$'+wlabel,fontsize=18)
	# ax.set_xlim(0,klim)
	# ax.set_ylim(0,wlim)
	im = plt.imshow(tr,**kwargs,extent=extent,cmap=cmap,clim=clim)
	del tr
	if (cbar):
		fig.colorbar(im)
	# fig.savefig(home_path+'k_w',bbox_inches='tight')
	return fig, ax


def plotting(fig,ax,name,default='.png'):
	# default save 
	if default == '.png':
		second = '.jpeg'
	else:
		second = '.png'

	try:
		fig.savefig(name+default,bbox_inches='tight')
		print('plotted using {} format'.format(default))
	except:
		fig.savefig(name+second,bbox_inches='tight')
		print('plotted using {} format'.format(second))
	plt.clf()
	return None

def plotDispersion(transmatrix, klimlow, klimup, wlimlow, wlimup, cbar=False, clim = (None,None),  labels=False): 
	tr = np.log10(transmatrix)[1:,1:] # "[1:,1:]" gets rid of the first row and coloumn
	# often want to chop to get better colour contrast. So pass in the extent. Usually this will be [0, klim, 0, wlim] where klim and wlim come from "plotting_vals"
	extent=[klimlow,klimup,wlimlow,wlimup]
	# experiment with clim, "(None, None)" means it will choose it for you, usually what you want. Sometimes I alter it to get better contrast
	if (labels):
		plt.xlabel('Wavenumber' + r' $[\omega_{c}/V_{A}]$',fontsize='15') #need to change when normalising 
		plt.ylabel('Frequency' + r' $[\omega_{c}]$',fontsize='15') #need to change when normalising 
	im = plt.imshow(tr,**kwargs,extent=extent,cmap='jet',clim=clim)#, clim = (-4.0,None))
	if (cbar): # off by default
		plt.colorbar()
	del tr
	return im
	
# =================================================

# Finds the limits of a batch size with its start and stop parameters
def batchlims(n,batch_size,index_list,remainder):
	batch_ini, batch_fin = n*batch_size+1, (n+1)*batch_size # assures no overlap between start and stop positions between iterations
	if n == 0: # should account for if it is the first one, it includes 0th file and hence needs to be changed
		batch_ini = 0
	if remainder: # accounts for remainder, needs to be passed new argument which will be True or False
		if remainder==1:
			print('Remainder Handled with 0th file.')
		elif remainder==2:
			print('Remainder x1...')
			batch_ini = index_list[-1]
			batch_fin = batch_ini
			print('Handled.')
		else:
			print('Remainder...')		
			batch_ini = index_list[-1]+2-n
			batch_fin = index_list[-1]
			print('Handled.')
	return batch_ini, batch_fin

# returns the batched fieldmatrix files for a range of field quantities. Also loads one of them if "load" parameter is true
def get_batch_fieldmatrix(index_list,quantities=['Magnetic_Field_Bz'],quantity='Magnetic_Field_Bz',load=True,para=False):
	batch_size, StartStop = BatchStartStop(index_list) # dealing with large number of files
	times = np.zeros(len(index_list))
	meanBz = 0
	for i in range(StartStop.shape[0]):
		(start,stop)=StartStop[i,:]
		print('batch #: {} , start : {}, stop : {}'.format(i,start,stop))
		tSarr = np.arange(start,stop+1,1)
		Sarr = np.zeros((len(tSarr),2),dtype='object')
		for j in range(len(tSarr)):
			Sarr[j,:] = [tSarr[j], quantities]
		pool=mp.Pool(mp.cpu_count())
		fm = np.array(pool.map_async(getFm,Sarr).get(99999)) # parallel loading of the individual files, not the batches
		pool.close()
		times = batch_getTimes(times,start,stop)
#		if i == 0: meanBz = np.mean(fm[0:10])
#		fm = fm-meanBz
		for q in range(len(quantities)):
			dumpfiles(fm[:,q,:],'batch_'+str(i)+'_fieldmatrix_'+quantities[q])
	dumpfiles(times,'times')
	if load:
		fieldmatrix = load_batch_fieldmatrix(index_list,quantity,para=para)
		return times, fieldmatrix
	else:
		return times, None

def loadFm(arr):
	ind, quantity = arr
	return read_pkl('batch_'+str(ind)+'_fieldmatrix_Magnetic_Field_Bz')


def getFm(arr):
	ind, quantities = arr
	vals = []
	try: 
		nfile = sdf.read(('%05d'%ind)+'.sdf')
	except:
		try: 
			nfile = sdf.read(('%04d'%ind)+'.sdf')
		except: 
			return None # if cant load at all, will kick up errors later
	for q in range(len(quantities)):
		vals.append(getQuantity1d(nfile, quantities[q]))
	return vals

# batch loads the fieldmatrix of a given quantity
def load_batch_fieldmatrix(index_list=[],quantity='Magnetic_Field_Bz',para=False):
	os.chdir(os.getcwd()) ## refresh directory in case new files have been made
	quantities = getFields()
	## find number of batch_files, if None then load normally
	batch_lst = np.array([i for i in os.listdir() if 'batch_' in i and 'fieldmatrix_'+quantity in i])	
	batch_lst = [i for i in range(len(batch_lst))]
	print('batch_lst :: ',batch_lst)
	if len(batch_lst) !=0 : # dumped in batches
		print('batch size :: {}'.format(int(len(index_list)/len(batch_lst))))
		batch_size = int(len(index_list)/len(batch_lst))
		print(batch_lst)
		if para:
			loadfunc = loadFm
			print('Parallel Load:')			
			tbatch_lst = np.zeros((len(batch_lst),2),dtype='object')
			for j in range(len(batch_lst)):
				tbatch_lst[j,:] = [batch_lst[j], quantity]
			pool=mp.Pool(mp.cpu_count())
			fieldmatrix=np.vstack(np.array(pool.map_async(loadfunc,tbatch_lst).get(99999)))
			pool.close()
		else:
			try: # see if quantity is batched
				print('Loop Load:')
				fieldmatrix = []
				for num in batch_lst: ## TODO ; make NxM matrix then assign
					fieldmatrix.append(read_pkl('batch_'+str(num)+'_fieldmatrix_'+quantity))
				fieldmatrix = np.vstack(np.array(fieldmatrix))
			except:
				print('## ERROR ## : Couldnt load batch {}'.format(num))
				_, fieldmatrix = get_batch_fieldmatrix(index_list,quantities,quantity,load=True)
	else: # not dumped in batches
		try:
			fieldmatrix = read_pkl('fieldmatrix_'+quantity)
		except:
			_, fieldmatrix = get_batch_fieldmatrix(index_list,quantities,quantity,load=True)
	print('Fieldmatrix fully loaded')
	return fieldmatrix

# returns the angle of the magnetic field wrt x and y
def getMagneticAngle(d0):
	# get the angle the magnetic field makes to the x and y directions
	# returns the angle (1d if in z-x plane or 2d if in z-x-y volume)
	# returns the angle in degrees
#	d0 = sdfread(0)
	try:
		Bx = getMeanquantity(d0, 'Magnetic_Field_Bx')
		By = getMeanquantity(d0, 'Magnetic_Field_By')
		Btot = getMeanField3D(d0, 'Magnetic_Field_B')
		phi_x, phi_y = np.arccos(Bx/Btot), np.arccos(By/Btot)
	except:
		phi_x = float(input('B0 angle to xhat [deg]::'))*(const.PI/180)
		phi_y = const.PI/2
	return phi_x, phi_y # will return 90 degrees for phi_y most of the time

# Plots the cold plasma dispersion for ionic species 1 and 2 (two maj or maj and min)
def coldplasmadispersion(file0,omegas,theta=None): 
	# Assumes one of the species is always electrons (harcoded)
	# angle between B & x-hat
	# returns:
			# k1, k2, k3 solutions (not-normalised)	
	if not theta: 
		theta, _ = getMagneticAngle(file0) # assuming B directed in z-x plane
	else:
		theta = theta*const.PI/180
	sin = np.sin(theta) ; cos = np.cos(theta)
	print(theta, sin, cos)
	# setup electron plasma and cyc freq
	wpe = getPlasmaFreq(file0,'Electrons')
	wce = getCyclotronFreq(file0,'Electrons',getChargeNum('Electrons'))
	wpf = [wpe]
	wcf = [wce]
	# loop over all ion species
	species_lst = getIonSpecies(file0)
	for species in species_lst:
		if species == '':
			wps = 0
			wcs = 0
		else:
			wcs = getCyclotronFreq(file0,species,getChargeNum(species))
			wps = getPlasmaFreq(file0,species)
		wpf.append(wps)
		wcf.append(wcs)
	# convert to numpy arrays
	wpf=np.array(wpf)
	wcf=np.array(wcf)
	# setup components	
	l = len(omegas)
	R=np.ones(l); P=np.ones(l); L=np.ones(l); S=np.zeros(l); D=np.zeros(l) 
	B=np.zeros(l); F=np.zeros(l); A=np.zeros(l); C=np.zeros(l) 
	tr=0; tl=0; tp=0
	# calculate components
	for i in range(len(wpf)):
		tr += (wpf[i]**2)/(omegas*(omegas + wcf[i]))
		tl += (wpf[i]**2)/(omegas*(omegas - wcf[i]))
		tp += (wpf[i]**2)/(omegas**2)
	R = R - tr
	L = L - tl
	P = P - tp

	S = 0.5*(R+L) ; D = 0.5*(R-L)
	C = P*R*L
	B = R*L*(sin**2) + P*S*(1.0 +cos**2)
	F = (((R*L - P*S)**2)*(sin**4) + 4.0*(P**2)*(D**2)*(cos**2))**0.5
	A = S*(sin**2) + P*(cos**2)
	n1 = np.zeros(l, dtype=complex) ; n2=np.zeros(l, dtype=complex); n3=np.zeros(l, dtype=complex)#; n4=np.zeros(l, dtype=complex) 
	n3 = np.lib.scimath.sqrt((R*L)/S)
	n1 =  np.lib.scimath.sqrt((B+F)/(2.0*A))
	n2 = np.lib.scimath.sqrt((B-F)/(2.0*A))
	n3 = -np.lib.scimath.sqrt((B+F)/(2.0*A))
	#n4 = -np.lib.scimath.sqrt((B-F)/(2.0*A))
	del R, P, L, S, D, B, F, A
	return (np.real(n1)*omegas)/const.c , (np.real(n2)*omegas)/const.c , np.real((n3*omegas)/const.c) #, (n4*omegas)/c, omegas
	
# Plots the power spectrum of a signal from a 2d FFT (trans) matrix
def powerspectrum(trans,wnorm,wklims=[None,None],wkmax=[0,None,0,None]):
	# Plots Fourier power as a function of frequency.
	# trans - input matrix that has been 2D FT'd
	# wklims - obtain from "getDispersion" limits, the limit (normalised) of the total FFT
	# wkmax = harmonicmin,harmonicmax,kmodelow,kmodehigh - maximum that the power is taken over (normalised)
	# if you want the power between the 1st and 10th harmonics, harmonicmin = 1,harmonicmax = 10 # normalisation dependent
	wlim,klim = wklims
	harmonicmin, harmonicmax, kmodelow, kmodehigh = wkmax	
	nw,nk = trans.shape
	wstart = int(nw*(harmonicmin/wlim)) 
	wstop = int(nw*(harmonicmax/wlim))
	print(wstart,wstop)
	## Wavenumber range
	kstop = int(nk*(kmodehigh/klim))
	kstart = int(nk*(kmodelow/klim))
	print(kstart,kstop)
	#kstop = int((trans.shape[1]/2.0) + (kmodehigh/klim)*(trans.shape[1]/2.0))
	#kstart = int((trans.shape[1]/2.0) - (kmodelow/klim)*(trans.shape[1]/2.0))
	power = np.zeros((wstop-wstart))
	omegas = np.zeros((wstop-wstart))
	dw = wnorm*wlim/nw # assuming all equally spaced
	for i in range(wstart,wstop):
		power[i-wstart] = np.sum((trans[i,kstart:kstop])**2) #positive k only
		omegas[i-wstart] = (i+wstart)*dw
#	omegas = (wlim*wnorm)*np.linspace((wstart/trans.shape[0]),(wstop/trans.shape[0]),len(power))
	return np.log10(power), omegas #notice it returns the log of the power and un-normalised (rad/s) omegas

# Same as the Fourier power obtained in powerspectrum() but this summates in k instead
def powerspectrum_k(trans,wlim,klim,harmonicmin,harmonicmax,kmodelow,kmodehigh):
	wstart = int(trans.shape[0]*(harmonicmin/wlim))
	wstop = int(trans.shape[0]*(harmonicmax/wlim))
	kstop = int(((kmodehigh/klim))*trans.shape[1])
	kstart = int(((kmodelow/klim))*trans.shape[1])

	power_k = np.zeros((kstop-kstart))
	
	for i in range(kstart,kstop):
		power_k[i-kstart] =  np.sum((trans[wstart:wstop,i])**2)
	
	wavenumbers = np.linspace((kstart*klim/trans.shape[1]),(kstop*klim/trans.shape[1]),len(power_k))

	return np.log10(power_k), wavenumbers


# Plots the change in field energy densities from their mean value
def getEnergies(energy_quant,fieldquant,nt,dump=True):
	F0 = np.zeros(len(energy_quant))
	for i in range(len(energy_quant)):
		if energy_quant[i] == 'Magnetic_Field_Bz':
			F0[i] = getQuantity1d(sdfread(0), 'Magnetic_Field_Bz')

	Energies = np.zeros((len(fieldquant),nt))
	Energies_mat = np.zeros((len(fieldquant),nt,len(getGrid(sdfread(0))[0]))) # energies, time, space
	for t in range(0,nt):
		if t%(nt//20)==0: print(str(round(100*t/nt))+'...%') # print every 5%
		d = sdfread(t)
		for s in range(len(fieldquant)):
			if 'Field' in fieldquant[s]:
				Energies_mat[s,t,:] = (getQuantity1d(d,fieldquant[s])-F0[s])**2
				Energies[s,t] = np.mean(Energies_mat[s,t,:])#(getQuantity1d(d,fieldquant[s])-F0[s])**2)
			else:
				Energies_mat[s,t,:], Energies[s,t] = getTotalKineticEnergyDen(d,fieldquant[s])

	## energy_quant and fieldquant should be the same length 
	if dump:
		for i in range(len(energy_quant)):
			dumpfiles(Energies[i,:],energy_quant[i])
			dumpfiles(Energies_mat[i],energy_quant[i]+'matrix')
			
	return Energies_mat, Energies


def getEnergyfromFieldmat(fieldmatrix,fieldquant):
	(nt,nx) = fieldmatrix.shape
	if fieldquant == 'Magnetic_Field_Bz':
		F0 = np.mean(fieldmatrix[0:10,:]) # mean Bz so have deltaBz
	
	Energy = np.zeros(nt)
	for t in range(0,nt):
		Energy[t] = np.mean((fieldmatrix[t,:]-F0)**2)

	return Energy

	
# Gets the total kinetic energy density to be used in plotenergies()
def getTotalKineticEnergyDen(d,species):
	#"Derived_EkBar" gives the KE per grid cell avergaed over all particles. In this function I then multiply by the density of particles to get the total
	#EK density. Ideally, would have KE for each particle, sum, then divide by "volume". Can calculate this from EPOCH momenta output 
	# (see function in "possibly_useful_funcs.py", needs modifying), but uses a lot of disc space and RAM hungry. Will probably have to do this for publication though...
	try:
		ek = getQuantity1d(d, "Derived_EkBar_"+species)
	except:
		ek = getQuantity1d(d, "Derived_Average_Particle_Energy_"+species)

	den = getQuantity1d(d, "Derived_Number_Density_"+species)
	energy_density = ek*den
	mean_energy_density = np.mean(energy_density)

	return energy_density, mean_energy_density

### Makes a video out of images with the image extension given
#def makevideo(image_foler=os.getcwd(),video_name='video',image_ext='.png',video_ext='.gif'):
#	images = [img for img in os.listdir(image_folder) if img.endswith(image_ext)]
#	frame = cv2.imread(os.path.join(image_folder, images[0]))
#	height, width, layers = frame.shape
#	video = cv2.VideoWriter(video_name+video_ext, 0, 1, (width,height))
#	for image in images:
#		video.write(cv2.imread(os.path.join(image_folder, image)))
#	cv2.destroyAllWindows()
#	video.release()
#	return None	

## Dump pkl files
def dumpfiles(array, quant):
	print('Pickling '+quant+'...')
	with open(quant+'.pkl', 'wb') as f:
		pickle.dump(array,f)

## Read pkl files
def read_pkl(quant):
	with open(quant+'.pkl', 'rb') as f:
		print('Loading '+quant+'...')
		array = pickle.load(f)
	print('Done.')
	# automatically closed when loaded due to "with" statement
	return array

## Scans dir for pkl files with "name" and returns boolean if it's in dir
def scan_pkl(names):
	# could implement a removal of "name" based off of returned boolean array: 
	# https://stackoverflow.com/questions/14537223/remove-items-from-a-list-using-a-boolean-array
	# names can be a single str or an array of files to load
	try: #check size of array
		np.shape(names)[0] 
		names=names
	except:
		names=[names]
	indir=[] # boolean array
	for name in names:
		try:
			with open(name+'.pkl', 'rb') as f: dummy = True # if can open then its in the dir
			indir.append(True)
		except:
			indir.append(False)
	if False in indir: #pkl file not in dir therefore need to run check again
		return False
	else: # all pkl files given in dir so can continue
		return True
		
# Calculating the LH freq using an effective mass of the majorities rather than the classical calculation
def LowerHybridMassEffective(file0,wpe,wce,n_e):
	m_eff = getEffectiveMass(file0)
	w_PI2 = (n_e*const.qe**2)/(m_eff*const.e0)
	w_CI2 = (const.qe*getMeanField3D(file0,'Magnetic_Field_B')/m_eff)**2
	W_LH = 1/np.sqrt((1/w_PI2)+1/(wce*np.sqrt(w_CI2)))
#	W_LH = (((w_PI2)+(w_CI2))/(1+((wpe**2)/(wce**2))))**.5	
	return W_LH # not normalised

def getUpperHybrid(wce,wpe):
	return np.sqrt(wce**2+wpe**2)

def getLowerHybrid(wpi,wci,wpe,wce):
	return np.sqrt(((wpi**2)+(wci**2))/(1+((wpe**2)/(wce**2))))

## Calculates the cold wave modes for a given plasma (LH, UH, Omode, Xmode, Light, Cyclotron harmonics, CA) # doesnt actually calculate FAW, this is done separately
def ColdWaveModes(ax,omega_lst,va,wnorm,OMODE=False,UH=False,LH=False,CA=False,CYC=False,LGHT=False,Lcut=False,Rcut=False,LHact=False,W1=False,W2=False):
	wpi, wpe, wci, wce = omega_lst
	w = wnorm*np.linspace(0,50,10000) # normalised in units of wci and va ## TODO; change this to be more general
	k = (wnorm/va)*np.linspace(0,200,10000)
	## O-mode
	if OMODE:
		wom = np.sqrt(wpe**2+(const.c**2)*(k**2))
		plt.plot(k*va/wnorm,wom/wnorm,color='k',linestyle='--')#,label='O-mode')
	## Upper-hybrid
	if UH:
		wuh = getUpperHybrid(wce,wpe)
		plt.axhline(wuh/wnorm,color='k',linestyle='-.')#,label='UH')
	## Lower-hybrid ; full description, typically use LHact
	if LH:
		wlh = getLowerHybrid(wpe,wci,wpe,wce)
		plt.axhline(wlh/wnorm,color='k',linestyle='--')#,label='LH') ## handled in the FAW calculation
	## Above is also LH but full version, this is the simplified version	
	if LHact:
		wlh = 1/np.sqrt((1/(wpi**2))+(1/(wce*wci)))
		plt.axhline(wlh/wnorm,color='k',linestyle='--')
	## Compressional Alfven
	if CA:
		wA = k*va
		plt.plot(k*va/wnorm,wA/wnorm,color='k',linestyle='--')#,alpha=0.5,label='Comp Alf')
	## Light travelling in a vacuum
	if LGHT:
		w = k*const.c
		plt.plot(k*va/wnorm,w/wnorm,linestyle='--',color='k')#label='light')
	## Cyclotron Harmonics
	if CYC:
		if wnorm == wci:
			for i in range(0,16):
				plt.axhline(i,color='k',linestyle=':',alpha=0.5)
		elif wnorm == wce:
			for i in range(0,4):
				plt.axhline(i,color='k',linestyle=':',alpha=0.5)
		else:
			None
	## Buchsbaum 1960 2 ion resonance, W1 & W2
	if W1:
		w1 = np.sqrt(wce**2 + wpe**2)
		plt.axhline(w1/wnorm,color='k',linestyle=':')
	if W2:
		w2 = np.sqrt(wce*wci*((wce*wci + wpe**2)/((wce**2)+(wpe**2))))
		plt.axhline(w2/wnorm,color='k',linestyle='-.')
	## L cut-off frequency
	if Lcut:
		sqt = np.sqrt((wci-wce)**2+4*(wpe**2+wpi**2+wce*wci))
		wLcut1 = 0.5*((wci-wce)+sqt)
		wLcut2 = 0.5*((wci-wce)-sqt)
		if wLcut1 > 0:
			plt.axhline(wLcut1/wnorm,color='white',linestyle='-.')
		if wLcut2 > 0:		
			plt.axhline(wLcut2/wnorm,color='white',linestyle=':')
	## R cut-off frequency
	if Rcut:
		sqt = np.sqrt((wci-wce)**2+4*(wpe**2+wpi**2+wce*wci))
		wRcut1 = 0.5*((wce-wci)+sqt)
		wRcut2 = 0.5*((wce-wci)-sqt)
		if wRcut1 > 0:
			plt.axhline(wRcut1/wnorm,color='white',linestyle='-.')
		if wRcut2 > 0:		
			plt.axhline(wRcut2/wnorm,color='white',linestyle=':')
	return ax

## Returns all normalisation values used later (frequencies and lengths)
def getFreq_Wavenum(file0,maj_species,min_species,Zmaj,Zmin):
	wci  = getCyclotronFreq(file0,maj_species,Z=Zmaj)
	wce  = getCyclotronFreq(file0,'Electrons',Z=1)
	va   = getAlfvenVel(file0)
	lambdaD= getDebyeLength(file0,'Electrons')
	wpe	 = getPlasmaFreq(file0,species='Electrons')
	wpi  = getPlasmaFreq(file0,species=maj_species)

	return wci,wce,va,lambdaD,wpe,wpi


def Chi0Calc(file0,v0,kall,omegaall,species='Deuterons',wci=None,theta=90):
#	Adapted from appendix A in ref. https://aip.scitation.org/doi/10.1063/1.860304 ; eqns [A.6], [A.8] and [A.16]
	# IN # 
	# v0 : perpendicular birth energy of minority species 
	# kall : real k2 solutions to the cold plasma dispersion
	# omegall : corresponding (to kall) omega solutions along the FAW 
	# species : species name you want to normalise frequencies to
	# theta : Angle B0 to sim domain
	# OUT #
	# Chi0 : Solution to eqn. [29]
	# PIxx, PIxy, PIyy : Components to the solution for Chi0, given in the ref		 

	if not wci:
		wci = getCyclotronFreq(file0,species)
	theta = theta*(const.PI/180) # radians

	PIxx = np.zeros(len(omegaall),complex) ; PIxy = np.zeros(len(omegaall),complex) ; PIyy = np.zeros(len(omegaall),complex)

	for i in range(0, omegaall.shape[0]):
		l = round(omegaall[i]/wci) #l closest to the omega
		k = kall[i]
		kpara = kall[i]*np.cos(theta)
		kperp = kall[i]*np.sin(theta)
		
		za = kperp*v0/wci
		zarr = np.linspace(0,2*za,10000)
		dzarr = (zarr[-1]-zarr[0])/len(zarr)

		# try and except clause to remove l=0 division
		try:
			PIxx[i] = spec.jv(2*l,2*za)
		except:
			None
		try:
			PIxy[i] = (-1j*za/l)*(spec.jvp(2*l,2*za))
		except:
			None
		try:
			PIyy[i] = (1-(za/l)**2)*spec.jv(2*l,2*za)+(za/(2*l**2))*np.sum(dzarr*spec.jv(2*l,zarr))
		except:
			None

	vA = getAlfvenVel(file0)
	Chi02 = PIxx - ((2*1j)*omegaall*wci/((kall**2)*(vA**2)))*PIxy + ((wci**2)/((vA**2)*(kall**2)))*PIyy
	Chi0 = np.sqrt(Chi02)
	
	return Chi0#, PIxx, PIxy, PIyy


## Find the growth rates of the MCI in its linear phase based off of drift and spread velocities
def growth_rate_man(minions, majions, theta, file0, u, vd, vr, karr, omegarr):
	# vr is para drift
	# vd is perp drift
	# emin is beam energy in eV

	# THIS GROWTH RATE CORRESPONDS TO EQUATIONS (8)-(10) OF MCCLEMENTS ET AL. POP 3 (2) 1996
	# I DON'T USE EQ. (11) TO GET THE FREQUENCY, RATHER I USE THE COLD PALSMA DISPERSION RELATION
	
	Zmin = getChargeNum(minions)
	Zmaj = getChargeNum(majions)
	wcycmin = getCyclotronFreq(file0,minions) # cyc freq for beam ions (minority)
	wcyci = getCyclotronFreq(file0,majions) # cyc freq for bulk ions

	wpmin = getPlasmaFreq(file0, minions) #square of beam ion plasma frequency, z=1 for tritons
	wpi = getPlasmaFreq(file0, majions) #square of bulk ion plasma frequency, z=1 for deuterons
	vA = getAlfvenVel(file0)

	theta = theta*(const.PI/180.0) #radians
	gammas = np.zeros(omegarr.shape[0]) #growth rates

	for i in range(0, omegarr.shape[0]):
		l = round(omegarr[i]/wcycmin) #l closest to the omega
		k = karr[i]
		kpara = karr[i]*np.cos(theta)
		kperp = karr[i]*np.sin(theta)

		Npara = (kpara*vA)/omegarr[i]
		Nperp = (kperp*vA)/omegarr[i]

		eetal = (omegarr[i] - kpara*vd - l*wcycmin)/(kpara*vr) # vd=0
		za = kperp*u/wcycmin
		
		w = omegarr[i]
		Jl = spec.jv(l,za)
		Jlprime = spec.jvp(l,za)
		JlJlprime = Jl*Jlprime

		############## N_l ###############################################################
		nlterm1 = -1*2*l*(w/wcyci)*(JlJlprime/za)
		nlterm2 = ((w**2-wcyci**2)/wcyci**2)*((Npara**2)*(((l**2 * Jl**2)/za**2)+Jlprime**2)+Nperp**2 * ((l**2 * Jl**2)/za**2))
		nlterm3 = ((l**2 * Jl**2)/za**2) + Jlprime**2
		nl = nlterm1 + nlterm2 + nlterm3
		##################################################################################
						#
		############# M_l #################################################################
		mlterm1 = 2*l*(w/wcyci)*(Jlprime**2 + (Jl**2/za**2)*(l**2 - za**2))
		mlterm2 = -1*2*((w**2-wcyci**2)/wcyci**2)*(JlJlprime/za)*(l**2 * Nperp**2 - (za**2 - 2*l**2)*Npara**2)
		mlterm3 = 2*(JlJlprime/za)*(za**2 - 2*l**2)
		ml = mlterm1 + mlterm2 + mlterm3
		#################################################################################		
						#
		########## Gamma ################################################################
		gterm1 = (wpmin**2/wpi**2)*(wcyci**4/((wcyci+(w-wcyci)*Npara**2)*(wcyci-(w+wcyci)*Npara**2)))
		gterm2 = ((l*wcycmin/(kpara*vr))*ml-((2*u**2)/(vr**2)*eetal*nl))
		gterm3 = np.sqrt(const.PI)/(2*w)*np.exp(-1*eetal**2)
		gammas[i] = gterm1*gterm2*gterm3
		#################################################################################		
	
	####### Want gamma > 1 only ###########
	posomega = [] ; posgamma = []	

	for i in range(0,gammas.shape[0]):
		if (gammas[i] >= 0):
			posomega.append(omegarr[i])
			posgamma.append(gammas[i])

	return np.array(posomega,dtype='float'), np.array(posgamma,dtype='float')
	
	
## Find and plot the distribution function of the species specified
def dist_fn(index_list,xyz,species):
	# Returns the distribution function plotted as a KDE from the binned data of the particles momenta, converted to velocity then normalised in units of the species thermal velocity
	# index_list :: same input as usual
	# xyz :: dimension 'x' , 'y' or 'z' depending on which component you want to measure
	# species :: name of the species you want to measure, typically take the minority ion
	# returns a 2d array of velocities through time and space for each file which contains the correct information, shape (nfiles, space)
	if xyz == 'x':
		mom = 'Particles_Px_'+species
		vel = r'$v_x$'
		vname = 'Vx_'
	elif xyz == 'y':
		mom = 'Particles_Py_'+species
		vel = r'$v_y$'
		vname = 'Vy_'
	elif xyz == 'z':
		mom = 'Particles_Pz_'+species
		vel = r'$v_z$'
		vname = 'Vz_'
	else:
		print('# ERROR # : Invalid velocity option.')
		raise SystemExit

	v = []
	vth = getThermalVelocity(0,species) # assumes 0th file has temperature
	mass = getMass(species)
	for i in index_list:
		d = sdfread(i)
		try:
			v.append((getQuantity1d(d,mom)/mass)/vth)
		except:
			continue
	del d

	nval = 1000
#	bins = 750
	fig,ax=plt.subplots(figsize=(8,5))
	print('Plotting KDEs...')
	for v_arr in v:
		print('vel file {}...'.format(i))
		kde = stats.gaussian_kde(v_arr)
		vv = np.linspace(min(v_arr),max(v_arr),nval) # nval is just the length of the plotted array in v and f(v)
		#	ax.hist(v_arr,bins=bins,density=True)    # change number of bins accordingly, density=True normalises the histogram
		ax.plot(vv,kde(vv))

	ax.set_xlabel(vel+r'$/v_{th}$',fontsize=18)
	ax.set_ylabel(r'$f(\mathbf{v})$',fontsize=18)
	fig.savefig(vname+species+'_dist_fn.jpeg',bbox_inches='tight') # plots all KDEs on the same axis
	
	dumpfiles(v,vname+species)
	return v

################# Bicoherence #################

def getBicoh(karea, warea, signal, dt, T , L, wnorm, knorm, nfft, noverlap, window=False, bispectrum=False,klabel=r'$\lambda_{De}$',cmap='jet'):
	
	print('ka :: {}\nwa :: {}'.format(karea,warea))
	verts = [
	   (0., 0.),  # left, bottom
	   (0.,karea),  # left, top
	   (warea, karea),  # right, top
	   (warea, 0.),  # right, bottom
	   (0., 0.),  # ignored
	]

	codes = [
		Path.MOVETO,
		Path.LINETO,
		Path.LINETO,
		Path.LINETO,
		Path.CLOSEPOLY,
	]
	
	area = Path(verts, codes)
#	try:
#		bicoh = read_pkl('Bicohmat_karea_{}_warea_{}'.format(karea,warea))
#		extent1 = [0, karea, 0, karea]
#	except:
	extent1 = [0, karea, 0, karea]
	if bispectrum:
		bispec,bicoh=bispectrum2D(signal, dt, L, T, area, nfft, noverlap, norm=[knorm,wnorm], window=window, bispectrum=bispectrum) 	
		dumpfiles(bispec,'Bispecmat_ka_wa_{}_{}_nfft_{}_noverlap_{}'.format(karea,warea,nfft,noverlap))
	else:
		bicoh=bispectrum2D(signal, dt, L, T, area, nfft, noverlap, norm=[knorm,wnorm], window=window, bispectrum=bispectrum)
	
	dumpfiles(bicoh,'Bicohmat_ka_wa_{}_{}_nfft_{}_noverlap_{}'.format(karea,warea,nfft,noverlap))

	if bispectrum: # written to do both bicoh and bispec at the same time
		fig, ax, _, _ = plot_bicoh(bispec, extent=extent1, bispectrum=True, smooth=True, cbar=True, cmap=cmap) ## bispec
		ax.set_xlabel(r'$k_1$'+klabel,fontsize=18)
		ax.set_ylabel(r'$k_2$'+klabel,fontsize=18)
		plotting(fig,ax,'bispectrum_karea_{}_warea_{}'.format(karea,warea))
		fig1, ax1, _, _ = plot_bicoh(bicoh, extent=extent1, bispectrum=False, smooth=True, cbar=True, clim=(0, 1), cmap=cmap) ## bicoh
		ax1.set_xlabel(r'$k_1$'+klabel,fontsize=18)
		ax1.set_ylabel(r'$k_2$'+klabel,fontsize=18)
		plotting(fig1,ax1,'bicoherence_karea_{}_warea_{}'.format(karea,warea))
	else: # dont know why you would use this when bispectrum=True does both
		fig, ax, _, _ = plot_bicoh(bicoh, extent=extent1, bispectrum=True, smooth=True, cbar=True, clim=(0, 1), cmap=cmap)
		ax.set_xlabel(r'$k_1$'+klabel,fontsize=18)
		ax.set_ylabel(r'$k_2$'+klabel,fontsize=18)
		plotting(fig,ax,'bicoh_karea_{}_warea_{}'.format(karea,warea))

	print('MEAN BICOH :: ',np.mean(bicoh))
	return fig, ax


# New cold dispersion for waves brought on by the Lower Hybrid Instability (LHDI)
def calcNewColdDisp(in_klimprime,Te,Ti):
	# Calc and Plot new dispersion relation with modified Lower-Hybrid (LH) freq
	# eqns (5) and (9) from this paper ; http://aip.scitation.org/doi/10.1063/1.3132628
	# Tobias S-H 05/10/2022 
	file0 = sdfread(0)
	ion_species, ion2_species, min_ions = getIonSpecies(file0)
	
	B = getMeanField3D(file0, 'Magnetic_Field_B') # T 
	me = const.me
	c = const.c
	mi = getMass(ion_species)
	m_min = getMass(min_ions)
#	Ti = 1E3 * eV_to_K # keV ==> K
#	Te = 1E3 * eV_to_K # keV ==> K
	print(Te, Ti)
	Eth_e = Te * const.kb # Joules
	theta, _ = getMagneticAngle(file0)
	theta = theta
	n0 = getMeanquantity(file0, 'Derived_Number_Density_Electrons')
	ni = getMeanquantity(file0, 'Derived_Number_Density_'+ion_species)
	n_min = getMeanquantity(file0, 'Derived_Number_Density_'+min_ions) 
	
	w_pi = getPlasmaFreq(file0,ion_species)
	w_pe = getPlasmaFreq(file0,'Electrons')
	w_ci = getCyclotronFreq(file0,ion_species)
	w_ce = getCyclotronFreq(file0,'Electrons')
	w_LH = 1/np.sqrt((1/(w_pi**2))+(1/(w_ci*w_ce)))
	vA = getAlfvenVel(file0)
	knorm = w_ci/vA
	Nk = 2500
	k = knorm*np.linspace(0,in_klimprime,Nk)
	Ve = np.sqrt(2*Eth_e/me) # electron thermal speed

	w_LH_new = np.sqrt((w_LH**2) * (1 + (np.cos(theta)**2)*(mi/me)))
#	plt.axhline(w_LH_new/w_ci, color='k', alpha=0.5, linestyle='--')

	W1 = 3*(Ti/Te)*(1+(w_pe**2)/((c*k)**2)) + (w_pe**2)/(2*(k*c)**2) + 9/2 - (15+21*(w_pe/k*c)**2)/(4*(1+(w_pe/k*c)**2)**2)
	W2 = (-3*((w_pe/(k*c))**2) + (1-6*(w_pe/k*c)**2)/(4*(1+(w_pe/k*c)**2)**2)) * (mi/me) * (np.cos(theta))**2
	W3 = 3*(((mi/me)*(np.cos(theta)**2) + ((w_pe/(k*c))**2) - (Ti/Te))/(1 + (w_pe/k*c)**2 + (mi/me)*np.cos(theta)**2))
	W = W1 + W2 + W3

	omega = np.sqrt((w_LH_new**2)/(1+(w_pe/(k*c))**2)*(1 + ((mi/me)*np.cos(theta)**2)/((1 + (w_pe/(k*c))**2)) + W*(k*Ve/w_ce)**2))
#	plt.plot(k/knorm,omega/w_ci, color='k', linestyle='-')
#	plt.ylabel('w/w_ci',fontsize=18)
#	plt.xlabel('kv_A/w_ci',fontsize=18)
#	plt.show()
	return k, omega, w_LH_new


def calcPowerICE(FT_2d, wlim, klim, wcyc, va, kwidth, wmax, kmax):
	# IN:
	# FT_2d : whole matrix of spatiotemporal FFT
	# wlim : normalised limit on freq
	# klim : normalised limit on wavenumber
	# wcyc : cyclotron freq used for normalisation
	# va : Alfven velocity
	# kwidth : width (normalised) to take power
	# wmax, kmax : maximum limits to take power up to
	# OUT:
	# power : summed 
	# omegas : array of frequencies from 0 to wmax and len FT_2d.shape[0]
	(nw,nk) = FT_2d.shape
	cellim = [nw*(wmax/wlim), nk*(kmax/klim)]
	FT_2d = FT_2d[:int(cellim[0]),:int(cellim[1])] # within max lims
	FT_2d = FT_2d[1:,1:] # remove first lines 
	print(nw,nk,cellim)
	kwidth = nk*kwidth/klim
#	plt.imshow(np.log10(FT_2d),aspect='auto',interpolation='nearest',cmap='magma',origin='lower')
	power = []
	kcell=[]
	for i in range(FT_2d.shape[0]): # range over w
		kcell.append(np.where(FT_2d[i,:]==max(FT_2d[i,:])))
		kcell[-1] = (kcell[-1])[0]
		obv = 5*(nk/klim)
		if i > 0: # can only check against previous kcells therefore require at least one point already in array
			if kcell[-1] > kcell[-2]+obv:
				FT_sorted = np.sort(FT_2d[i,:])
				kcell[-1] = np.where(FT_2d[i,:]==FT_sorted[-2])
				kcell[-1] = (kcell[-1])[0] # find second largest value
#		plt.scatter(kcell[-1],i,color='k',s=10)
		krange_min = kcell[-1] - kwidth
		krange_max = kcell[-1] + kwidth
		print(krange_max, krange_min)
		if krange_min < 0:
			krange_min = 0
		if krange_max > nk*kmax/klim:
			krange_max = cellim[1]-2
#		plt.scatter(krange_max,i,color='red',s=10)	
#		plt.scatter(krange_min,i,color='blue',s=10)
		power.append(np.sum(FT_2d[i,int(krange_min):int(krange_max)]**2))

	omegas=np.linspace(0,wmax,FT_2d.shape[0])
	plt.plot(omegas,power)
	plt.yscale('log')
	plt.show()
	
	return power, omegas

#def getPoynting(times,nt,nx,min_species='Alphas',plot=True):
#	sval=[]
#	for i in range(0,nt): # find restart files with all Electric and Magnetic Fields
#		try:
#			_ = getQuantity1d(sdfread(i),'Electric_Field_Ey') # see if restart file or not
#			sval.append(i)
#		except:
#			None
#	print(sval)
#	# vector components
#	Sx = np.zeros((len(sval),nx))
#	Sy = np.zeros((len(sval),nx))
#	Sz = np.zeros((len(sval),nx))
#	for j in range(len(sval)):	
#		print(str(np.around(100*j/len(sval),2))+'%')
#		d = sdfread(int(sval[j]))
#		Ex = getQuantity1d(d,'Electric_Field_Ex')
#		Ey = getQuantity1d(d,'Electric_Field_Ey')
#		Ez = getQuantity1d(d,'Electric_Field_Ez')
#		Bx = getQuantity1d(d,'Magnetic_Field_Bx')
#		By = getQuantity1d(d,'Magnetic_Field_By')
#		Bz = getQuantity1d(d,'Magnetic_Field_Bz')
#		Sx[j,:] = (1/const.mu0)*(Ey*Bz - Ez*By)
#		Sy[j,:] = -1*(1/const.mu0)*(Ex*Bz-Ez*Bx)
#		Sz[j,:] = (1/const.mu0)*(Ex*By - Ey*Bx)
#	# remove first line
#	Sx = Sx[1:,:]
#	Sy = Sy[1:,:]
#	Sz = Sz[1:,:]
#	Smag = np.sqrt(Sx*Sx + Sy*Sy + Sz*Sz) # magnitude of the Poynting Vector
#	if plot:
#		plotPoynting(times,Smag,'mag',min_species=min_species) # plot magnitude of Poynting Vector by default
#	return Sx, Sy, Sz, Smag

#def plotPoynting(times,Smat,dim,cmap='magma',cbar=True,min_species='Alphas'):
#	tc_a = 2*const.PI/getCyclotronFreq(sdfread(0),min_species)
#	L = getGridlen(sdfread(0))
#	fig,ax=plt.subplots(figsize=(6,6))
#	im = ax.imshow(np.log10(Smat),**kwargs,cmap=cmap,extent=[0,L,0,times[-1]/tc_a])
#	ax.set_xlabel(r'$x$'+'  '+'[m]',fontsize=18)
#	ax.set_ylabel(r'$t/\tau_{c\alpha}$',fontsize=18)
#	if cbar: fig.colorbar(im, label=r'$log_{10}(\mathbf{S}_{mag})$')

#	if dim == 'x':
#		plt.savefig('Sx.png')
#		plt.clf()
#	elif dim == 'y':
#		plt.savefig('Sy.png')
#		plt.clf()
#	elif dim == 'z':
#		plt.savefig('Sz.png')
#		plt.clf()
#	elif dim == 'mag':
#		plt.savefig('Smag.png')
#		plt.clf()
#	else:
#		print('Not a suitable dimension for Poynting Vector.')
#		raise SystemExit


#def Poynting2dFT(times,nt,nx,L=None,in_klimprime=100,plot=True):
#	try:
#		if not L:
#			L = getGridlen(sdfread(0))
#		_,_,_,Smag = getPoynting(times,nt,nx,min_species='Alphas',plot=False)
#		FT2d = get2dTransform(Smag)
#		T = times[-1]
#		if plot:
#			d0 = sdfread(0)
#			klim = 0.5*2*const.PI*nx*1/L
#			wlim = 0.5*2*const.PI*Smag.shape[0]/T # modified freq limit
#			wnorm = getCyclotronFreq(d0,'Deuterons')
#			vA = getAlfvenVel(d0)
#			knorm = wnorm/vA
#			fig,ax = plt.subplots(figsize=(8,2))
#			ax.imshow(np.log10(FT2d),**kwargs,extent=[0,klim/knorm,0,wlim/wnorm])
#			ax.set_xlim(0,in_klimprime)
#			ax.set_ylim(0,wlim/wnorm)
#			ax.set_xlabel(r'$kv_A/\Omega_D$',fontsize=18)
#			ax.set_ylabel(r'$\omega/\Omega_D$',fontsize=18)
#			plotting(fig,ax,'FT_2d_Smag.png')
#	#		plt.show()
#		return FT2d
#	except:
#		return None


def Endict(field):
	Energydict = {
	'Magnetic_Field_Bz': 'Bzenergy',
	'Magnetic_Field_By': 'Byenergy',
	'Magnetic_Field_Bx': 'Bxenergy',
	'Electric_Field_Ez': 'Ezenergy',
	'Electric_Field_Ey': 'Eyenergy',
	'Electric_Field_Ex': 'Exenergy',
	}
	EnergyLabeldict = {
	'Magnetic_Field_Bz': r'$\langle\Delta B_z^2\rangle/2\mu_0$',
	'Magnetic_Field_By': r'$\langle B_y^2\rangle/2\mu_0$',
	'Magnetic_Field_Bx': r'$\langle B_x^2\rangle/2\mu_0$',
	'Electric_Field_Ez': r'$\langle E_z^2\rangle\epsilon_0/2$',
	'Electric_Field_Ey': r'$\langle E_y^2\rangle\epsilon_0/2$',
	'Electric_Field_Ex': r'$\langle E_x^2\rangle\epsilon_0/2$',
	}
	Energymult = {
	'Magnetic_Field_Bz': const.mult_magn_field,
	'Magnetic_Field_By': const.mult_magn_field,
	'Magnetic_Field_Bx': const.mult_magn_field,
	'Electric_Field_Ez': const.mult_elec_field,
	'Electric_Field_Ey': const.mult_elec_field,
	'Electric_Field_Ex': const.mult_elec_field,
	}

	return Energydict.get(field), EnergyLabeldict.get(field), Energymult.get(field)

## get energy quantities, names, labels and multipliers
def getEnergyLabels(file0,species):
	## check which fields are in normal dump
	fieldquant = getFields()
	## append to energy_quant array
	energy_quant=[] ; names=[] ; Energy_mult=[]
	for field in fieldquant:
		fieldlabel, labelval, Emult = Endict(field)
		energy_quant.append(fieldlabel)
		names.append(labelval)
		Energy_mult.append(Emult)
	
	## append ion arrays
	ionquant = [] ; ionlabels = []
	for spec in species:
		if spec != '':
			ionquant.append(spec+'_KEdens')
			ionlabels.append(getIonlabel(spec))
			fieldquant.append(spec)
	
	## append all 
	for k in range(len(ionquant)):
		energy_quant.append(ionquant[k])
		names.append(ionlabels[k])
		Energy_mult.append(1) # no multiplier for particles

	return energy_quant, names, Energy_mult, fieldquant
	
def getFields(n=1):
	## check if field values in normal dump 
	try:
		d = sdfread(n) # hardcoded for now, will default if cant load # TODO
		keys = getKeys(d)
		fieldquant = []
		for key in keys:
			if 'Field' in key:
				fieldquant.append(key) # append if is a field value
	except:
		fieldquant=['Magnetic_Field_Bz','Magnetic_Field_By','Electric_Field_Ex'] # default fields
	
	return fieldquant

# Home function which can be called, will call on getEnergies(). Generates KE_Dens and field E_Dens for species 1, 2, 3 and e-
def energies(sim_loc,frac=1,plot=False,leg=True,integ=False,linfit=False,electrons=True):
	width = 10
	height = 6
	#height = width*(1/const.g_ratio)
	fig,ax=plt.subplots(figsize=(width,height))
	
	d0 = sdfread(0)
	species = list(getIonSpecies(d0))

	if species[-1] == '': # check if has minority ion
		species[-1] = getAllSpecies(d0)[-1] # check if minority exists as negative ion

	maj_species, maj2_species, min_species = species
	if electrons: species = np.append(species,'Electrons')

	energy_quant, names, Energy_mult, fieldquant = getEnergyLabels(d0,species)
	
	## do when no minority ring-beam
	if maj2_species == '': 
		if min_species == '': # Single ion
			Single = True
			Double = Triple = False
			print('Single Ion...')
		else: # 2 ions:
			Double = True
			Single = Triple = False
			print('Two Ions...')
	else: # 3 ions
		Triple = True
		Single = Double = False
		print('Three Ions...')

	Zmaj1 = getChargeNum(maj_species)
	Zmaj2 = getChargeNum(maj2_species)
	Zmin = getChargeNum(min_species)

	index_list = list_sdf(sim_loc)
	os.chdir(sim_loc)

	# loop :: species & fields
	# parallel loading and dump total
	try:
		times=read_pkl('times')
	except:
		times = batch_getTimes(np.zeros(len(read_pkl('Exenergy'))),1,len(read_pkl('Exenergy'))-1) # janky approach, only used for testing and when times.pkl hasn't been made
	n = (len(index_list)//frac)
	
	#################################
	csum = 0

	# Calculate
	print(energy_quant)
	if set([i+'.pkl' for i in energy_quant]).issubset(set(os.listdir(os.getcwd()))):
		print('All energy -pkl- files present.')		
	else:
		energy_data_mat, energy_data = getEnergies(energy_quant,fieldquant,n,dump=True) # mat,_ = getEnergies()
#		for i in range(energy_data.shape[0]):
#			dumpfiles(energy_data[i,:],energy_quant[i])

	if plot:
		print('Plotting energies...')
		colors=['b','cyan','g','r','m','orange','k','salmon','lightgreen'] # will only use all of them if there are 3 +ve and 1 -ve species ## assuming no extra field values
		tnorm=2*const.PI/getCyclotronFreq(d0,min_species) # last species
		
		mean_to = 10
		dt = (times[-1]-times[0])/len(times)
		print('### dt :: ',dt)

		left, bottom, width, height = [0.5,0.15,0.3,0.22]
#		ax2 = fig.add_axes([left,bottom,width,height]) # inset of zoomed in portion of energy dens
		for i in range(0,len(energy_quant)): # loops over all species given in energy_mult
			Energy=read_pkl(energy_quant[i])
			mean_Energy=np.mean(Energy[:mean_to])
			energy_plot = (Energy-mean_Energy)*Energy_mult[i]
			ax.plot(times[::frac]/tnorm,energy_plot[::frac],label=names[i],color=colors[i])
#			ax2.plot(times[::frac]/tnorm,energy_plot[::frac],color=colors[i])
#			ax2.set_xlim(0,0.1)
			if integ:
				csum += integrate(energy_plot,dt) # integral of each species, sums over all species
			if linfit:
				thresh = times/tnorm < 1
				A = np.vstack([times[thresh]/tnorm, np.ones(len(times[thresh]/tnorm))]).T
				m, c = np.linalg.lstsq(A, energy_plot[thresh], rcond=None)[0]
				#	m = np.mean(np.gradient(Energy[thresh],times[thresh]/tnorm))
				print(m)
				Energy = Energy-(m*times/tnorm+c)
				ax.plot(times/tnorm,Energy,label=labels[i]+'  '+' m='+str(m//1),color=colors[i])
	
		print('### FINAL INTEGRAL :: {}'.format(csum))
		timeLabel = r'$t$'+getOmegaLabel(min_species)+r'$/2\pi$'
		ax.set_xlabel(timeLabel,fontsize=18)
		ax.set_ylabel(r'$\Delta u$'+'  '+'['+r'$Jm^{-3}$'+']',fontsize=18)
		if (leg): 
			ax.legend(ncol=1,loc='best',fontsize=14,labelspacing=0.1,borderpad=0.1) # change ncol to make legend span multiple columns
		ax.set_xlim(0,max(times/tnorm))
		ax.axhline(0,color='darkgrey',linestyle='--')
#		ax.set_xlim(0,0.1)
		plotting(fig,ax,'energy_densities') # check to see if can save in .jpeg or .png (Orac is old)
		print('Plot saved.')
		plt.clf()
		del times, Energy, energy_plot
	
	return csum

def integrate(data,dt):
	csum = 0
	for dval in data:
		csum += dval*dt
	return csum

# Home function for calculating power spectrum (call when in the relevant dir)
def power(wnorm,wklims=[None,None],wkmax=[None,None],norm_omega=r'$\Omega_D$',quantity='Magnetic_Field_Bz',plot=False,read=True,outp=True):
	wlim_prime,klim_prime=wklims
	wmax,kmax=wkmax
	if read:
		try:
			log10_power = read_pkl('log10_power')
			omegas = read_pkl('omegas_power')
			calc = False
		except:
			print('Can\'t read power and omegas...')
			read = False ; calc = True
	else:
		calc = True
	if calc:
		try:
			FT_2d = read_pkl('FT_2d_'+quantity)
		except:
			print('## ERROR ## :: Cannot load FT_2d')
			raise SystemExit
		print('Calculating power...')
		log10_power,omegas=powerspectrum(FT_2d,wnorm,[wlim_prime,klim_prime],[0,wmax,0,kmax])
		dumpfiles(log10_power,'log10_power')
		dumpfiles(omegas,'omegas_power')
	
	if plot:
		print('Plotting Power...')
		width,height = 8,5
		fig,ax=plt.subplots(figsize=(width,height))

		ax.plot(omegas/wnorm,10**log10_power)
		for i in range(1,int(wmax)+1):
			ax.axvline(i,color='k',alpha=0.5,linestyle=':')
		ax.set_xlabel(r'$\omega$'+'/'+norm_omega,fontsize=18)
		
		ax.set_xlim(0,wmax)
		ax.set_ylabel('Power',fontsize=18) ## unitless
		## find min, max
		paxmin, paxmax = min(10**log10_power[10:]), max(10**log10_power[10:])
		ax.set_ylim(paxmin/10,10*paxmax) ## power of 10 higher and lower than min / max
		ax.set_yscale('log')
		plotting(fig,ax,'power_'+quantity)

	if outp: 
		return omegas, log10_power
	else: 
		del omegas; del log10_power ; return None


def BatchStartStop(ind_lst,default=700):
	if len(ind_lst) < default: # choose
		batch_size = len(ind_lst)
	else:
		lower = 100 ; upper = default + 101		
		batch_size = [i for i in range(lower,upper,1) if (len(ind_lst)-1)%i==0] # create array of working batch sizes
		if (len(batch_size)) == 0: # no batches available
			batch_size = default # assign a default value
		else:
			batch_size = batch_size[-1] # choose largest applicable
	print('batch_size :: {}'.format(batch_size))
	start=[i*batch_size for i in range(0,len(ind_lst)) if i*batch_size<len(ind_lst)] # get starting positions of batches
	StartStop = [[i,i+batch_size-1] for i in start if i+batch_size-1 < len(ind_lst)] # get start and end positions of batches
	if StartStop[-1][1] % ind_lst[-1] != 0: # remainder
		if (ind_lst[-1] - StartStop[-1][0]) > 1024: # above file limit, need to split again
			StartStop.append([StartStop[-1][1]+1,ind_lst[-1]]) # very rarely this is the case
		else:# below file limit, can assign a larger final batch size
			StartStop[-1][1] = int(ind_lst[-1]) # typically just adds one to the end
	return batch_size, np.array(StartStop)

## concentration ratios, number of species dependent
def getConcentrationRatios(file0):
	species = getIonSpecies(file0)
	n0 = getMeanquantity(file0,'Derived_Number_Density_Electrons')
	xiarr = []
	for spec in species:
		try:
			xiarr.append(getMeanquantity(file0,'Derived_Number_Density_'+spec)/n0)
		except:
			xiarr.append(0)
	return np.array(xiarr) # xi1, xi2, xi3

## returns the signed curvature of a set of datapoints
def signedCurv(fx,dx):
	### take the absolute for defined "curvature" or the square integral to get the total curvature of fx 
	## In
	#	fx : np array of points to calculate curvature
	#	dx : spacing in x-array points
	## Out
	#	kx : signed curvature, as defined here (https://en.wikipedia.org/wiki/Curvature#Graph_of_a_function)
	#	kxint : Integral of squared signed-curvature, higher the value --> more curvy / less smooth
	fxp = np.gradient(fx,dx)
	fxpp = np.gradient(fxp,dx)
	kx = abs(fxpp)/(1+fxp**2)**(1.5) # reduce to just "curvature"
	kxint = np.sum(kx*dx)
	return kx, kxint
	
def paraVelocity(INDEX):
	index, species, mSpec = INDEX
	return np.sqrt(2*getQuantity1d(sdfread(index),'Derived_Average_Particle_Energy_'+species)/mSpec)


# quote "cigarette plots", which show the distribution of velocity/energy normalised to the  
def ciggies(sim_loc,species_lst=['Deuterons','Alphas'],nval=10000,para=False,eload=True,vload=False,logo=False):
	if logo: # print the ciggies logo # made this because I could
		import pyfiglet
		print(pyfiglet.figlet_format('===#  --  @@@',font='letters',justify='center',width=110))
		del pyfiglet

	if eload:
		vload = False
	elif vload:
		eload = False
	elif eload == vload:
		print('## ERROR ## :: vload and eload are set equal') ; raise SystemExit
	else:
		print('## ERROR ## :: no value set for eload vload') ; raise SystemExit

	# initial parameters, dimensions and constants
	time = read_pkl('times')
	d0 = sdfread(0)
	Nx = len(getGrid(d0)[0])
	Nt = len(time)
	vA = getAlfvenVel(d0)
	species_lst = [i for i in species_lst if i != ''] # remove missing species
	
	# setup mass and velocity arrays
	for species in species_lst:
		L = getGridlen(d0)
		tcSpecies = 2*const.PI/getCyclotronFreq(d0,species)
		T = time[-1]/tcSpecies
		dens = getMeanquantity(sdfread(0),'Derived_Number_Density_'+species)
		# try loading f first rather than instantly making it
		if vload: # velocities
			ylabel = r'$v/v_A$'
			cbar_label = r'$\log_{10}[f(v)]$'
			figname = 'fv_vA_'
			matnorm = vA
		if eload: # energies
			ylabel = r'$E$' + '  ' + '['+r'$keV$'+']'
			cbar_label = r'$\log_{10}[f(E)]$'
			figname = 'fE_keV_'
			matnorm = 1000*const.qe
			try:
				fSpecies = read_pkl(figname+species)
				fload = True
			except:
				fload = False
		# check if fload is possible
		if fload:
			if vload: # should load if here, can't be here otherwise
				matSpecies = read_pkl('v_'+species)
			if eload:
				matSpecies = read_pkl(species+'_KEdensmatrix')/dens
			# normalise
			matSpecies = matSpecies/matnorm
			matmin  = np.min(matSpecies) ; matmax = np.max(matSpecies) ; matmean = np.mean(matSpecies)
		else: # calculate distribution if can't load it already
			ind_lst = list_sdf(sim_loc)
			massSpecies = getMass(species)
			try: # try loading
				if vload:
					matSpecies = read_pkl('v_'+species)
				if eload:
					matSpecies = read_pkl(species+'_KEdensmatrix')/dens
			except: # didn't load, calculate instead
				if eload:
					matSpecies,_ = getEnergies([species+'_KEdens'],[species],Nt,dump=True)/dens # loads energy density matrix of species
				if vload:
					tind_lst = np.zeros((len(ind_lst),3),dtype='object')
					for i in range(len(ind_lst)):
						tind_lst[i,:] = [ind_lst[i], species, massSpecies]
					if para: # parallel calculation, just for velocities ## TODO: change this to be more general
						pool = mp.Pool(mp.cpu_count()//2)
						matSpecies = np.vstack(np.array(pool.map_async(paraVelocity,tind_lst).get(99999)))
						pool.close()
					else:
						for t in range(Nt):
							if int(100*t/Nt)%5==0:print(100*t/Nt, ' %')
							matSpecies[t,:] = np.sqrt(2*getQuantity1d(sdfread(t),'Derived_Average_Particle_Energy_'+species)/massSpecies)		
					dumpfiles(matSpecies,'v_'+species)
			# normalise
			matSpecies = matSpecies/matnorm
			matmin  = np.min(matSpecies) ; matmax = np.max(matSpecies) ; matmean = np.mean(matSpecies)
			# initialise distribution matrix
			fSpecies = np.zeros((Nt,nval))
			# calculate fv/fE matrix
			for t in range(fSpecies.shape[0]): # loop through time
				xarr = np.linspace(0,L,Nx)
				yarr = matSpecies[t,:]
				fSpecies[t,:],_,_ = np.histogram2d(xarr,yarr,range=[[0,L],[matmin,matmax]],bins=(1,nval),density=True)
			# dump distribution matrix
			dumpfiles(fSpecies,figname+species)

		# extents of matrix per species
		extents = [0,T,matmin,matmax]
		print('minmat ',matmin,'maxmat ',matmax)

		# setup figure
		fig,ax = plt.subplots(figsize=(8,8/const.g_ratio))
		im = ax.imshow(np.log10(fSpecies.T),**kwargs,extent=extents,cmap='jet')

		# colorbar and labels
		cbar = plt.colorbar(im)
		cbar.ax.set_ylabel(cbar_label, rotation=90,fontsize=18)
		ax.set_xlabel(r'$t$'+getOmegaLabel(species)+r'$/2\pi$',fontsize=16)
		ax.set_ylabel(ylabel,fontsize=16)
		
		# set y-lim
		if (0.9*matmean) < matmin:
			if (1.1*matmean) > matmax:
				ax.set_ylim(0.9*matmean,1.1*matmean)
			else:
				ax.set_ylim(0.9*matmean,maxmat)
		else:
			if (1.1*matmean) > matmax:
				ax.set_ylim(matmin,1.1*matmean)
			else:
				ax.set_ylim(matmin,matmax)

		plt.gca().ticklabel_format(useOffset=False)
		# savefig
		fig.savefig(figname+species+'.png',bbox_inches='tight')
	return None

# extract peaks in a dataset (used for power spectra comparisons)	
def extractPeaks(data,Nperwcd=1,prominence=0.3):
	return signal.find_peaks(data,distance=Nperwcd,prominence=prominence)[0] # tune till Nperwcd encapsulates all peaks (visually)	

# calculate the growth of k-modes for various time spacings
def map_k_growth(sim_loc,normspecies='Deuterons',omega_min=0,omega_max=100,domega=0.25,dt_frac=0.5,tstart_frac=0.5,tend_frac=2.0,color='k',plot=False):
	## define sim
	sim_loc = getSimulation(sim_loc)
	d0 = sdfread(0)
	## normalisation
	vA = getAlfvenVel(d0)	
	wnorm = getCyclotronFreq(d0,normspecies)
	klim = 0.5*2*const.PI/getdxyz(d0)
	times = read_pkl('times') ; tnorm = 2*const.PI/wnorm
	T = times[-1]/tnorm
	
	## cold plasma disp
	omegas = wnorm*np.arange(int(omega_min),int(omega_max),domega)
	species = getIonSpecies(d0)	
	#k1,k2,k3 = coldplasmadispersion(d0, omegas, theta=89) 
	k2 = omegas/(vA*np.sqrt(1+np.cos(86.3*const.PI/180.)**2))
	thresh = k2 > 0
	k2 = k2[thresh] ; omegas = omegas[thresh]
	## FT1d & extract data at karr points
	FT1d = read_pkl('FT_1d_Magnetic_Field_Bz')
	harmonics = FT1d.shape[1]*(k2/klim)
	# loop through harmonics, then loop through time window to acquire multiple growth rates for a given k-mode
	[t0,t1] = [tstart_frac*tnorm,tend_frac*tnorm]
	t1_2 = (t1+t0)/2
	dt = dt_frac*tnorm
	nmin = 1
	nmax = (t1+t0)/(2*dt) ; narr = np.arange(nmin,nmax+1,1)
	growthRates = np.zeros((len(harmonics),len(narr)))	
	for i in range(len(harmonics)):
#		data = np.mean((FT1d[1:,round(harmonics[i])-1:round(harmonics[i])+1]),axis=1)
		data = (FT1d[1:,round(harmonics[i])])
		# extract peaks in data
		peaks = extractPeaks(data,Nperwcd=10)
		datapeak = data[peaks]
		time_dummy = times[peaks]
#		indexpeak = np.ones(len(datapeak))*i
#		axs.scatter(indexpeak,time_dummy/tnorm,datapeak)
		# fit data over multiple times
		for n in range(len(narr)):
			try:
				thresh = (time_dummy > t1_2-narr[n]*dt) & (time_dummy < t1_2+narr[n]*dt)
				timefit = time_dummy[thresh]
				datafit = datapeak[thresh]
				popt, pcov = curve_fit(lambda t,b,c: np.exp(b*t)+c,timefit,datafit,maxfev=5000) # exponential fitting
				b,c = popt
				fitted = np.exp(b*timefit) + c
				growthRates[i,n] = 2*const.PI*b # not normalised
			except:
				continue
	growthRatesMean = np.mean(growthRates,axis=1)
	growthRatesSTD = np.std(growthRates,axis=1)
	if plot:
		thresh = growthRatesMean > 0
		growthRatesMean = growthRatesMean[thresh]
		omegas = omegas[thresh]
		plt.plot(omegas/wnorm,growthRatesMean/wnorm,'o-',color=color)
		plt.savefig('growth_kmodes.png')
	return omegas, growthRatesMean, growthRatesSTD

## Calculates the growth rates from the gradients of the energy densities for a quantity
def grad_energydens(simloc,normspecies='Deuterons',quant='Magnetic_Field_Bz',mean_to=10,N=300,conv2=True):
	times = read_pkl('times')
	dt = (times[-1]-times[0])/len(times)
	d0 = sdfread(0)
	wcD = getCyclotronFreq(d0,normspecies)
	tcD = 2*const.PI/wcD
	t_end = times[-1]
	if quant in getFields(): # field
		pklname,label,mult = Endict(quant)
	else: # particle 
		pklname = quant+'_KEdens'
		label = getIonlabel(quant)
		mult = 1
	Energyfield = read_pkl(pklname)
	meanEnergyfield = np.mean(Energyfield[:mean_to])
	Energyfield = (Energyfield-meanEnergyfield)*mult
	## running mean
	if conv2: # double-convolve
		Energyfield = np.convolve(Energyfield,np.ones(N)/N,mode='valid')
	# gradient of log value to get growth rates of exp growth [since d/dx exp(ax) = a exp(ax)] 
	gradu = np.convolve(np.gradient(np.log10(Energyfield),dt),np.ones(N)/N,mode='valid')
	timesplot = np.linspace(0,max(times),len(gradu))
#	plt.plot(timesplot/tcD,gradu*2*const.PI/wcD)
#	plt.ylabel(r'$\gamma/\Omega_D$',fontsize=20)
#	plt.xlabel(r'$t/\tau_{cD}$',fontsize=20)
	return timesplot, gradu


## Calculates the phase correlation between two signals, sig & sig0 (only need fft)
# returns the value of the shift corresponding to the phase difference between both arrays
def phaseCorrelation(sig,fft_sig0,dw,wnorm,wmax=35):
	fft_sig = np.fft.fft(sig)
	fft_sig_conj = np.conj(fft_sig)
	R = (fft_sig0 * fft_sig_conj) / abs(fft_sig0 * fft_sig_conj)
	r = np.fft.ifft(R)
	shift = dw*np.argmax(r)
	if shift > dw*len(r)/2:
		shift = shift-wmax*wnorm
	return shift


# Plots and shows the experiment vs theory plot for the ratio between two species change in energy density
# Also has the ability to plot du per-particle ratio (per species) through time for each simulation (nrows) -- will need to un-comment these lines
def majIons_edens_ratio(sims,species=['Deuterons','Tritons'],time_norm=r'$\tau_{cD}$',\
						xlabel=r'$[\Delta u_1/\Delta u_2]_{max}$',ylabel=r'$(\xi_1/\xi_2)(m_2/m_1)(q_1/q_2)^2$',labels=[1,11,50],lims=((0,1),(0,1))):
	mean_to = 10
	N=50
	c=0
#	fig, ax = plt.subplots(figsize=(10,len(sims)*3-1),nrows=len(sims),sharex=True)
#	fig.subplots_adjust(hspace=0.15)#hspace=0.
	plt.plot(lims[0],lims[0],color='darkgray',linestyle='--') # 1:1 line
	home = os.getcwd()
	for sim in sims:
		sim_loc = getSimulation(sim)
		times = read_pkl('times')
		d0 = sdfread(0)
		nx = len(getQuantity1d(d0,'Derived_Number_Density'))
		n0 = getQuantity1d(d0,'Derived_Number_Density_Electrons')
		wc1 = getCyclotronFreq(d0,species[0])
		tc1 = 2*const.PI/wc1
		dt = (times[-1]-times[0])/len(times)
		dt_p = dt/tc1
		## field & species energies
		maxdu = np.zeros(len(species))
		duarr=[]
		xi1,xi2,_=getConcentrationRatios(d0)
		xi2_xi1 = xi2/xi1
		print('Secondary conc.:: ',xi2)
		marr = [getMass(species[0]),getMass(species[1])]
		qarr = [getChargeNum(species[0]),getChargeNum(species[1])]
		for i in range(len(species)):
			Energypart = read_pkl(species[i]+'_KEdens')/(1000*const.qe) # keV/m^3
			meanEnergypart = np.mean(Energypart[:mean_to])
			Energypart = np.convolve(Energypart,np.ones(N)/N,mode='valid')
			timespart = np.linspace(0,max(times),len(Energypart))
			thresh = timespart/tc1 < 6.1
			timespart = timespart[thresh]
			du = (Energypart[thresh]-meanEnergypart)
			maxdu[i] = np.max(du)
			duarr.append(du)
		plt.scatter((xi1/xi2)*(marr[0]/marr[1])*((qarr[1]/qarr[0])**2),maxdu[1]/maxdu[0])
#		duarr = np.array(duarr)
#		ax[c].plot(timespart/tcD,np.abs((xi2_xi1)*duarr[1]/duarr[0]))
#		ax[c].annotate(str(labels[c])+'%',(times[-1]/tc1-2,1e2),xycoords='data',**tnrfont)
#		ax[c].set_ylim(1e-3,1e3)
#		ax[c].axhline(marr[0]/marr[1],linestyle='--',color='k')
#		ax[c].set_yscale('log')
		os.chdir(home)
		c+=1
#	ax[0].set_yticks([1e-2,1e-1,1e0,1e1,1e2,1e3])
#	ax[1].set_yticks([1e-2,1e-1,1e0,1e1,1e2,1e3])
#	ax[2].set_yticks([1e-3,1e-2,1e-1,1e0,1e1,1e2,1e3])
#	midax = len(ax)//2
#	ax[0].set_xlim(0,times[-1]/tc1-1)
#	ax[int(midax)].set_ylabel(r'$\left(\frac{\xi_1}{\xi_2}\right)\left|\frac{\Delta u_1(t)}{\Delta u_2(t)}\right|$',fontsize=24)
#	ax[-1].set_xlabel(r'Time,  '+time_norm,**tnrfont)
#	fig.savefig('/storage/space2/phrmsf/dump/du_ratio_vs_time_primary.png')
	plt.ylabel(xlabel,**tnrfont) ; plt.xlabel(ylabel,**tnrfont)
	plt.xlim(lims[0]) ; plt.ylim(lims[1])
	plt.savefig('du_peak_vs_theory.png',bbox_inches='tight')
	plt.show()


def shared_area(sig1,sig2,fitgauss=False):
	lx = len(sig1)
	if lx!=len(sig2): 
		print('## ERROR ## :: Both signals need the same length')
		raise SystemExit
	dx = sig1[-1]/lx
	rolled = np.zeros(lx)
	sharea=[]
	for i in range(lx):
		rsig1 = np.roll(sig1,i)
		c1 = rsig1 < sig2
		c2 = sig2 < rsig1
		rolled[c1] = rsig1[c1]
		rolled[c2] = sig2[c2]
		sharea.append(np.sum(rolled*dx))
	if fitgauss:
		x = np.linspace(0,lx,lx)
		popt, pcov = curve_fit(lambda sharea,b,c: np.exp(a*(x-b)**2)+c,x,sharea,maxfev=5000) # exponential fitting
		peak = popt[1] # a,b,c 
		return sharea, peak
	return sharea
	
def outside_ticks(fig):
	for i, ax in enumerate(fig.axes):
		ax.tick_params(axis='both',direction='out',top=False,right=False,left=True,bottom=True)

def boutside_ticks(lax):
	for ax in lax:
		ax.tick_params(axis='both',direction='out',top=False,right=False,left=True,bottom=True)
	
def xoutside_ticks(lax):
	for ax in lax:
		ax.tick_params(axis='x',direction='out',top=False,right=False,left=False,bottom=True)

def ignorex(lax):
    for ax in lax:
        ax.tick_params(labelbottom=False)

def ignorey(lax):
    for ax in lax:
        ax.tick_params(labelleft=False)

## Sobel and Scharr kernel on an image which will return the gradient array
def Kernel(img,kernel='scharr'):
	# convert img to 0-255 color
	img = 255*img/np.nanmax(img)
	if kernel == 'scharr':
		Gx=np.array([[3,0,-3],[10,0,-10],[3,0,-3]])
		Gy=np.array([[3,10,3],[0,0,0],[-3,-10,-3]])
	if kernel == 'sobel':
		Gx=np.array([[1,0,-1],[2,0,-2],[1,0,-1]])
		Gy=np.array([[1,2,1],[0,0,0],[-1,-2,-1]])
	# magnitude and angle arrays
	kGmag = np.zeros(img.shape)
	kGangle = np.zeros(img.shape)
	# position relative to centre cell
	for i in range(img.shape[0]-1):
		for j in range(img.shape[1]-1):
			NORTH = SOUTH = WEST = EAST = True
			# four conditions on square image
			if i+1 > img.shape[0]: # reached end of img EAST
				EAST = False
			if i-1 < 0: # reached end of img WEST
				WEST = False
			if j+1 > img.shape[1]: # reached end of img NORTH
				NORTH = False
			if j-1 < 0: # reached end of img SOUTH
				SOUTH = False
			##
			if NORTH: N = img[i,j+1]
			else: N = 0
			#
			if SOUTH: S = img[i,j-1]
			else: S = 0
			#
			if WEST: W = img[i-1,j]
			else: W = 0
			#
			if EAST: E = img[i+1,j]
			else: E = 0
			#
			if NORTH and WEST: NW = img[i-1,j+1]
			else: NW = 0
			#
			if NORTH and EAST: NE = img[i+1,j+1]
			else: NE = 0
			#
			if SOUTH and WEST: SW = img[i-1,j-1]
			else: SW = 0
			#
			if SOUTH and EAST: SE = img[i+1,j-1]
			else: SE = 0
			## 
			timg = np.array([[NW,N,NE],[W,img[i,j],E],[SW,S,SE]])
			kG = np.array([np.sum(Gx*timg),np.sum(Gy*timg)])
			kGmag[i,j] = np.sqrt(kG[0]**2 + kG[1]**2)
			kGangle[i,j]= np.arctan2(kG[1],kG[0]) # radians between +- pi
	kGangle = np.nan_to_num(kGangle,posinf=np.nan,neginf=np.nan) # change +-inf vals to nan
	kGmag = np.nan_to_num(kGmag,posinf=np.nan,neginf=np.nan) 	 # change +-inf vals to nan
	return kGmag, kGangle
	
def checkallFields(ind_lst,quantities,quant):
	allFields = True
	for i in range(0,ind_lst[-1]):
		if set(getFields(i)) != set(quantities):
				allFields*=False
		else:
				allFields*=True
	allFields = bool(allFields)
	if not allFields: quantities=[quant]
	return quantities

def getSpectrogram(fieldmatrix,times,majspec='Deuterons',minspec='Alphas',nfo=[10,2],cellfreq=1,clim=(None,None),plot=True,cbar=True):
## in 
	# fieldmatrix :
	# times : 
	# majspec : 
	# minspec : 
	# nfo : list of integers with which to calculate nfft and noverlap (nfo[0] = len(times)/nfft, nfo[1]=nfft /noverlap)
	# cellfreq : 
	d0	= sdfread(0)
	tnorm = 2*const.PI/getCyclotronFreq(d0,minspec)
	nt,nx = fieldmatrix.shape

	# sampling freq
	fs = int(len(times)/(times[-1]-times[0]))
	# number of ffts to take and their overlap
	nfft = len(times)//nfo[0] # lower value --> better time resolution (worse freq res)
	noverlap = nfft//nfo[1] # better windowing between frames, removes bleeding
	# initialise matrices
	nes = (nt-nfft)//noverlap
	Sarr = []
	Spower = np.zeros((nfft,nx))
	# loop through each position and append Spower
	for i in range(0,nx,cellfreq): # loop through
		fmx = fieldmatrix[:,i]
		freqs, time, Sxx = signal.spectrogram(fmx,fs,nperseg=nfft,noverlap=noverlap)
		Sarr.append(Sxx)
		for j in range(len(time)):
			Spower[j,i] = np.sum(Sxx[:,j]**2)
	Sarr = np.array(Sarr)
	dumpfiles([Sarr,freqs,time],'Sxftarray_nfft_{}_noverlap_{}'.format(nfft,noverlap))
	dumpfiles(Spower,'Sxpower_nfft_{}_noverlap_{}'.format(nfft,noverlap))	

	if plot:	plotSpectrogram(Spower,times,tnorm,nfft,noverlap,minspec,clim,cbar)
	return Spower, Sarr, freqs, time

def plotSpectrogram(Spower,times,tnorm,nfft,noverlap,minspec='Alphas',clim=(None,None),cbar=True):
	Spower = Spower[np.all(Spower!=0,axis=1)]
	fig,ax=plt.subplots(figsize=(8,6))
	im = ax.imshow(Spower,extent=[0,1,0,times[-1]/tnorm],**kwargs,clim=clim)
	if cbar:
		cbar = plt.colorbar(im)
		cbar.set_label('Spectrogram Power',**tnrfont)
	ax.set_xlabel(r'$x/L$',**tnrfont)
	ax.set_ylabel(r'$t$'+getOmegaLabel(minspec)+r'$/2\pi$',**tnrfont)
	fig.savefig('Sx_mat_nfft_{}_noverlap_{}.png'.format(nfft,noverlap),bbox_inches='tight')
	return None