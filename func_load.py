import os, sys
# un-comment if PYTHONPATH variable isn't appended
path = os.getcwd()+'/lib'
sys.path.append(path)

# def f():
# 	global list_new
# 	global bispectral_analysis
# 	global my_constants
# 	global numpy
# 	global sdf
# 	global pickle
# 	global matplotlib
# 	global scipy
# 	global multiprocessing

from lib.list_new import *
from lib.bispectral_analysis import *
import lib.my_constants as const
# other packages
import numpy as np
import sdf
import pickle
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
#plt.style.use('classic')
#plt.tight_layout()
kwargs, tnrfont = formatting()

# confirmation message
print('Core modules imported.')

from scipy import special as spec
from scipy import stats, signal
from scipy import interpolate
from scipy.optimize import curve_fit
from matplotlib.path import Path
import multiprocessing as mp # for parallelisation
# confirmation message
print('Additional modules imported.')



'''
	Should be able now to make arrays using np.XXXX, plot using plt.XXXX and load 
	functions _without_ calling list_new.func .
	This needs to be called to each script used for analysing EPOCH.
'''
