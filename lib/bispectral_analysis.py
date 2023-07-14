#from functions import *
import numpy as np
from scipy.ndimage.filters import gaussian_filter
import matplotlib.pyplot as plt
import sys

#### NOTE: ALL THESE ROUTINES CALCULATE AND WORK WITH b^2, NOT b. ####

pi = np.pi

def calcNyquist2D(length, nx, duration, nt, angular=True):
	"""
	INPUT
	length: The length of the signal in the spatial dimension.
	nx: Number of data points in the spatial domain.
	duration: The duration of the signal.
	nt: Number of data points in the time domain.
	angular: True if one wants to output in rads/s, False if one wants 1/s.

	OUTPUT
	knyq: Nyquist frequency for space
	fnyq: Nyquist frequency for time
	"""

	deltaf = 1.0/duration
	if (angular):
		deltaf = 2.0*pi*deltaf

	deltak = 2.0*pi/length
	knyq = deltak*0.5*float(nx)
	fnyq = deltaf*0.5*(float(nt)-1)

	return deltak, deltaf, knyq, fnyq


def find_nearest_index(array,val):
	idx = (np.abs(array-val)).argmin()
	return idx, array[idx]
	

def bispectrum2D(data, dt, length, duration, area, nfft, noverlap=None, norm=[1,1], window=False, bispectrum=False):
	"""
	INPUT
	data: The 2D data matrix of shape (time, space) whos bicoherence or bispectrum we are to calculate.
	dt: The time between successive data points along the 0th dim of "data".
	length: The distance between successive data points along the 1st dim of "data".
	area: A PATH object which contains the segment of frequency and wavenumber space that we're interested in.
		For instance if the area is a square with k limits 0 to kmax, and frequency limits 0 to wmax, this routine will 
		search for couplings between all modes contained within this square. See this webpage for PATH documentation.
		https://matplotlib.org/api/path_api.html
	nfft: The width of the successive FT's.
	noverlap: Optional. The overlap of the successive FT's. Defaults to nfft//2.
	norm: Optional. norm[0] is for normalisation for wavevectors, norm[1] is for normalisation for frequencies.
	      Defaults to no normalisation .
	window: Optional. True uses a Hanning window in time, False uses no window.
	bispectrum: Optional. By default this routine calculates the bicoherence, but can return the bispectrum also.

	OUTPUT
	Ouputs the bicoherence only if bispectrum=False. If bispectrum=True, outputs bispectrum followed by bicoherence

	COMMENTS
	This routine is RAM hungry, do yourself a favour and chop/stride through your "data" matrix as much as you can before passing it in,
	as well as thinking carefully about "area".
	"""
	if (nfft > data.shape[0]):
		raise ValueError('not enough data points to perform an FFT of size %5i, aborting.' %nfft)
	if noverlap is None:
		noverlap = nfft // 2
	elif noverlap >= nfft:
		raise ValueError('noverlap must be less than nfft.')

	if bispectrum:
		print("Calculating bispectrum.")
	else:
		print("Calculating bicoherence.")

	nt = data.shape[0]
	nx = data.shape[1]
	deltak, deltaf, knyq, fnyq = calcNyquist2D(length, nx, duration, nt) # calculates nyq values for whole grid, need nyq values for size of FFT 

	deltak = deltak/norm[0]
	knyq = knyq/norm[0]
	deltaf = (2.0*pi/(dt*nfft))/norm[1]
	fnyq = 0.5*deltaf*nfft	

	klist = np.linspace(0, knyq, nx//2)
	#kmin = area.get_extents().ymin
	#kminind = find_nearest_index(klist, kmin)[0]
	kmax = area.get_extents().ymax
	kmaxind = find_nearest_index(klist, kmax)[0]
	nxinterest = 2*kmaxind # the maximum k of interest will be twice our largest k	
	
	nes = ((nt-nfft) // noverlap) + 1 #the number of FTs we'll do, might not cover all time
	newfnyq = (float(nfft)/float(nt))*fnyq
	
	if (window): #can make the window more flexible in the future
		window = np.outer(np.hanning(nfft), np.ones(nx))
	else:
		window = np.outer(np.ones(nfft), np.ones(nx))
	
	w, k = np.meshgrid(np.arange(-nfft/2, nfft/2)*deltaf, np.arange(-nxinterest/2, nxinterest/2)*deltak)  
	#for each k there are many w, viceversa. covering all 4 quadranrts of the Dispersion.
	w, k = -w.flatten(), k.flatten()
	#stack "nnft" lots of the same k on top of each other. stack "nx" lots of the same w on top of each other and reverse. 
	# length is nx*nfft.
	
	points = np.vstack((w, k)).T #shape is [nx*nfft,2], [:,0] are the above w's, [:,1] are the above k's.
	grid = area.contains_points(points) 
	# "area" is a path. "contains_points" returns whether the path contains points. boolean matrix, flat.
	grid = grid.reshape((nxinterest, nfft)).T  
	# The shape of this is nfft by nx. We now have the points in the w, k mesh that are in the enclosed "area".
	# Plotting this grid using imshow gives a polygon which "masks" the region of k,w space we're interested in.
#	plt.imshow(grid, origin='lower', aspect='auto', extent=[-fnyq, fnyq, -newfnyq, newfnyq]) ; plt.show() ; sys.exit(-1)
	
	kone = np.arange(nxinterest)
	ones = np.ones(nxinterest, dtype=np.int)
	kthree = np.zeros([nxinterest, nxinterest], dtype=np.int) - nxinterest # need to alter if want k1-k2. All elements are set to "-nx"
	kthree = kthree + np.outer(kone, ones) + np.outer(ones, kone) #k3=k1+k2 (3rd term is like k2, need to alter if want k1-k2
	# imagine the axis represent k1 and k2 (cells), and the value at a particular k1 and k2 represents k3 (cells). i.e kthree[n,m] = n+m
#	plt.imshow(kthree, extent=[-nx/2,nx/2,-nx/2,nx/2], origin='lower') ; plt.show() ; sys.exit()
	wthree = np.zeros([nxinterest,nxinterest], dtype=np.int) - nfft
	trans_interest_max = np.zeros(nxinterest, dtype='complex')

	numerator = np.zeros((nxinterest, nxinterest), dtype='complex')
	denominator = np.zeros((nxinterest, nxinterest))
	del w, k, points, kone	

	print("About to start loop. I have "+str(nes)+" iterations to do.")
## working version producing bicoh of wavenumbers k1 & k2
	for i in range(0, nes):
		print("Iteration: ", i+1)
		start = i*noverlap
		stop  = nfft + i*noverlap
		data_portion = data[start:stop,:] 

		trans=np.fft.fftshift(np.fft.fft2(data_portion*window))
		#plt.imshow(np.log10(np.abs(trans)), origin='lower', aspect='auto', extent =[-knyq, knyq, -fnyq, fnyq]) ; plt.show() ; sys.exit(-1)
		
		mid = trans.shape[1]//2
		trans = trans[:,mid-nxinterest//2:mid+nxinterest//2]
		trans_interest = grid*trans
#		plt.imshow(np.log10(np.abs(trans_interest)), origin='lower', aspect='auto') ; plt.show()# ; sys.exit(-1)
		trans_interest = np.fft.fftshift(trans_interest)
		trans = np.fft.fftshift(trans)
		# masking shifted transform with the area of interest, then shifting back (the mask was already shifted when created)

		wmax = np.argmax(np.abs(trans_interest), axis=0) #index of the value of omega for which this k has its max power
		for j in range(0,nxinterest):
			trans_interest_max[j] = trans_interest[wmax[j],j]
			# for each k, store the maximum power, based on what we just found.		

		wthree = np.outer(wmax, ones) + np.outer(ones, wmax) - nfft
		# w3=w1+w2 (3rd term is like w2, need to alter if want w1-w2

		bicoh = np.outer(trans_interest_max, trans_interest_max)*np.conj(trans[wthree,kthree]) 
		# F(f1)*F(f2)*[F(f3)*]

		numerator = numerator + bicoh
		denominator = denominator + np.abs(bicoh)
		del bicoh
	print("Finished loop")
	del wthree, trans_interest_max, wmax, trans_interest, data_portion, trans

	if (bispectrum): # half is empty, 1st element is Nan
		return (np.abs(numerator)[1:nxinterest//2, 1:nxinterest//2]), ((np.abs(numerator)/denominator)[1:nxinterest//2, 1:nxinterest//2])
	else:
		return (np.abs(numerator)/denominator)[1:nxinterest//2, 1:nxinterest//2]


def plot_bicoh(bicoh, extent=None, bispectrum=False, smooth=False, cbar=False, clim=(None, None), cmap ='jet'):
	if (smooth):
		bicoh = gaussian_filter(bicoh, sigma=1) 
	bicoh = np.transpose(np.tril(np.ones(bicoh.shape)))*bicoh
	bicohplot = np.ma.masked_where(bicoh < 0.001, bicoh) # 	bicohplot = np.ma.masked_where(bicoh < 0.001, bicoh)
#	cmap = plt.cm.jet # colormap used
#	cmap.set_bad(color='white') # sets a "bad" color, which corresponds to the mask
	if bispectrum:
		bicohplot = np.log10(bicohplot)
	fig,ax = plt.subplots()
	im = plt.imshow(bicohplot, interpolation='nearest', origin='lower', aspect='auto', extent=extent, cmap=cmap, clim=clim)
	if (cbar):
		cbar = plt.colorbar()
#		cbar.ax.tick_params(labelsize=14)

	return fig, ax, bicohplot, bicoh

def sumbicoh(bicoh, extent, cutoff=0, smooth=False):
	"""
	INPUT
	bicoh: The 2D bicoherence matrix which goes from 0 to kmax in both directions.
	       This can be b or b^2.
	extent: Extent of bicoherence matrix
	cutoff: Values below this wont be included in the sum. There is a lot of nonlinear "noise" in PIC sims.
	        Default is 0, i.e. everything is summed.
	smooth: Can choose to smooth the bicoh matrix in much the same way as is often done for plotting the bicoherence.
	        In general this results in a smoother b^2 plot.
	        ALWAYS check smoothing doesn't chanmge your conclusions. Default is no smoothing.

	OUTPUT
	sumb: 1D array containing the sum of the bicoherence (or its square). Dimensions 0 to 2*kmax.
	kvec: 1D array of wavenumbers from 0 to 2*kmax.
	"""

	if (smooth):
		bicoh = gaussian_filter(bicoh, sigma=1) 
	bicoh = np.transpose(np.tril(np.ones(bicoh.shape)))*bicoh # reorganise matrix as in the "plot_bicoh" function
	kvec = np.linspace(extent[0], 2.0*extent[1], 2*bicoh.shape[0]) #interactions up to kmax + kmax
	sumb = np.zeros(kvec.shape)
	for i in range(0, bicoh.shape[0]):
		for j in range(0, bicoh.shape[1]):
			b = bicoh[i,j]
			if (b > cutoff): # don't include "nonlinear noise"
				sumb[i+j] = sumb[i+j] + b
	return sumb, kvec


def load_bicohdat(fdat, fhead):
	"""
	INPUT
	fdat: Simple data field containing the 2D bicoherence/bispectrum data matrix
	fhead: Simple data file containing: nfft, noverlap, stridex, extent[0], extent[1], extent[2], extent[3]
		 where:
			 nfft: The width of the successive FT's
			 noverlap: The overlap of the successive FT's. Defaults to nfft//2
			 stridex: Number of points we skip in the simulation domain
			 extent: Corresponds to the bounds set by "area". Should already be normalised.
	
	OUTPUT
	b: bicoherence/bispectrum matrix
	extent: Extent of bicoherence/bispectrum matrix

	COMMENTS
	The structure of fhead is currently very simple. Change to a dictionary in the future.
	"""
	b = np.loadtxt(fdat, delimiter=',')
	header = np.loadtxt(fhead, delimiter=',')
	extent = [header[3], header[4], header[5], header[6]]
	return b, extent
	
