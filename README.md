
# README

This is the code used throughout the thesis of [Tobias Slade-Harajda](https://warwick.ac.uk/fac/sci/physics/research/cfsa/people/slade-harajda/) from 2021-2025. 

Various files are included to demonstrate their use in analysis of power spectra, energy densities, Fourier transformations, correlation and many more.

General use of the code will involve the list_new and batch_load programs which are found under /lib. These work to:
- `list_new` :: list of _all_ functions
- `batch_load` :: creates a simulation object which is then analysed and stored. Majority of first-run analysis is done here, then pre-processing is managed after files are made and stored

## Packages
We utilise standard python packages throughout, but make sure you have `pickle`. Pickle is a method of storing data which transfers large data-sets into manageable files which can be easily "dumped" and "loaded" using the `dumpfiles()` and `read_pkl()` functions found in `list_new`. Pickle is extremely useful due to it's flexibilty in storing data as it can even pkl images and re-load them. Note: this hasn't been used extensively in this thesis so please read their [documentation](https://docs.python.org/3/library/pickle.html).

Included also here are useful files from other Git users such as:
- Particles :: python particles used in plasma simulations, small scale but works like a gyro-kinetic code
- Polycoherence :: polycoherence used to calculate coherence between multiple signals in freq-space. Can be compared to bispectral analysis

There are various files used here that are old and will likely be removed so don't rely on their remaining here. Expecially those in the old/ sub-dir.

Other acknolwedgements will likely include Omstavan Samant, Ben Chapman and James Cook as the writing of a few of these functions was done by them and is treated as legacy code. Every single function has either been written (as is the case with the majority) or edited by myself to suit my needs better.


