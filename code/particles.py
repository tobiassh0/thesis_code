import time
import my_constants as const
import numpy as np
from plasmapy import formulary, particles
#from plasmapy.plasma import *
from plasmapy.particles import *
from plasmapy.particles.particle_class import valid_categories
from plasmapy.simulation.particletracker import ParticleTracker
from plasmapy.simulation import ParticleTracker
from plasmapy.plasma import Plasma
from plasmapy.formulary import ExB_drift, gyrofrequency
from astropy import units as u

deuteron = Particle('D+')
triton = Particle('T+')
alpha = Particle('alpha')
electron = Particle('e-')
print(electron.symbol,electron.mass)
particles = ParticleList([electron,deuteron,triton,alpha])

B0 = 2 * u.T
E0 = 0.5 * u.V/u.m
n_e = 1e19
T_e = 1000*const.qe/const.kb
v_perp = 1*u.m/u.s # np.sqrt(2*(T_e*const.kb)/const.me) * u.m/u.s
wCe = formulary.frequencies.gyrofrequency(B0, electron)
r_L = formulary.lengths.gyroradius(B0, electron, Vperp=v_perp)
print(wCe,r_L)

#start = time.time()
xarr=np.linspace(0,1,100) * u.m
yarr=xarr ; zarr=xarr # cube volume
#plasma = Plasma(
#	domain_x=xarr,
#	domain_y=yarr,
#	domain_z=zarr,
#)
#plasma.magnetic_field[2,...]=B0 # z-hat
#plasma.electric_field[0,...]=E0 # x-hat
##ExB_drift(plasma.electric_field[:,0,0,0],plasma.magnetic_field[:,0,0,0])

#t_end = 3*(2*np.pi*u.rad/wCe)
#print(t_end)
#nt = 100
#dt = t_end/nt
#trajectory = ParticleTracker(plasma, 'e-', 1, 1, dt, nt)
#trajectory.v[0][0] = v_perp
#trajectory.run()
#stop=time.time()
#print('Execution time :: {}'.format(str(stop-start)))
#trajectory.plot_trajectories()


#################################################
plasma = Plasma(
    domain_x=xarr,
    domain_y=yarr,
    domain_z=zarr,
)
B0 = 2 * u.T
plasma.magnetic_field[0, ...] = B0

E0 = 0.5 * u.V / u.m
plasma.electric_field[1, ...] = E0

ExB_drift(plasma.electric_field[:, 0, 0, 0], plasma.magnetic_field[:, 0, 0, 0])

freq = gyrofrequency(B0, "e-").to(u.Hz, equivalencies=u.dimensionless_angles())
gyroperiod = (1 / freq).to(u.s)
steps_to_gyroperiod = 50
timestep = 20 * gyroperiod / steps_to_gyroperiod

number_steps = steps_to_gyroperiod * int(2 * np.pi)
trajectory = ParticleTracker(plasma, "e", 1, 1, timestep, number_steps)

trajectory.v[0][0] = 1 * (u.m / u.s)
trajectory.run()
trajectory.plot_trajectories()


