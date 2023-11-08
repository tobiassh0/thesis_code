
# README

This is the code used throughout the thesis of [Tobias Slade-Harajda](https://warwick.ac.uk/fac/sci/physics/research/cfsa/people/slade-harajda/) from 2021-2025. 

Various files are included to demonstrate their use in analysis of power spectra, energy densities, Fourier transformations, correlation and many more.

General use of the code will involve the list_new and batch_load programs which are found under /lib. These work to:
- `list_new` :: list of _all_ functions
- `batch_load` :: creates a simulation object which is then analysed and stored. Majority of first-run analysis is done here, then pre-processing is managed after files are made and stored

## Packages
We utilise standard python packages throughout, but make sure you have `pickle`. Pickle is a method of storing data which transfers large data-sets into manageable files which can be easily "dumped" and "loaded" using the `dumpfiles()` and `read_pkl()` functions found in `list_new`. Pickle is extremely useful due to it's flexibilty in storing data as it can even pkl images and re-load them. Note: this has been used extensively in this thesis for fieldmatrix and Fourier transform storage, so please read their [documentation](https://docs.python.org/3/library/pickle.html).

Included also here are useful files from other Git users such as:
- [particles](https://github.com/PlasmaPy/PlasmaPy) :: taken from the PlasmaPy library. Simulates python particles used in plasma simulations, small scale but works like a gyro-kinetic code.
- [polycoherence](https://github.com/trichter/polycoherence) :: polycoherence used to calculate coherence between multiple signals in freq-space. Can be compared to bispectral analysis

There are various files used here that are old and will likely be removed so don't rely on their remaining here. Expecially those in the old/ sub-dir.

Acknolwedgements include Omstavan Samant, Ben Chapman and James Cook as the writing of a few of these functions was done by them and is treated as legacy code. Every single function has either been written (as is the case with the majority) or edited by myself to suit my needs better.


## Modular system
Have included the ability for the user to create a new file from scratch in the parent directory (i.e. above `/lib`). An example of this is listed below for a user who wants to analyse the spatial cross-correlation between multiple field components


```
# import all functions
from func_load import *
# import spatial correlation functions
import correlation.spatial_crosscor as sc
# instanstiate simulation location
sim = getSimulation(SIM_FILE_PATH)
# load array of times
times = read_pkl('times')
# list all available fields
fields = getFields()
# loop through all fields without repeating duplicates
for i in range(len(fields)):
    for j in range(i,len(fields)):
        F1 = fields[i]
        F2 = fields[j]
        _,_,norm1=Endict(F1)
        _,_,norm2=Endict(F2)
        dphase,Rtdx = sc.getSpatialCorrelation(d0=sdfread(0),F1=F1,F2=F2,norm=[norm1,norm2],plot=False)
        sc.plotSpatialCorrelation(F1,F2,Rtdx,dphase,times,MINSPECIES)        
```
