
from func_load import *

sim = getSimulation('/storage/space2/phrmsf/lowres_D_He3/0_10_p_90')
tcmin = 2*const.PI/getCyclotronFreq(sdfread(0),'Protons')
restart_files = np.arange(0,12500,500) # TODO; switch back to parallel checker
times=read_pkl('times')
print(restart_files)
species=['Deuterons','He3']
colors=['b','r']

# loop through species
for j in range(len(species)):
    try: # already dumped
        Vperp_squared = read_pkl('Vperpsquared_'+species[j])
        restart_times = np.array([times[t] for t in restart_files])
    except: # calculate
        print('Failed load Vperpsquared, calculating')
        # get Vperp and times at restart
        Vperp_squared,restart_times = linear_check_getVelocityPerp(restart_files,times,species[j],theta,theta_y)
        
    # RMS Vperp
    Vperp = np.sqrt(np.mean(Vperp_squared,axis=1)) # mean Vperp per time step
    Vperp_squared = Vperp**2

    mass_spec = getMass(species[j])
    tmagnetic_field = load_batch_fieldmatrix([],'Magnetic_Field_Bz')
    magnetic_field = np.zeros((len(restart_files),tmagnetic_field.shape[1]))
    for t in range(len(restart_files)):
        tindex = restart_files[t]
        magnetic_field[t,:] = tmagnetic_field[tindex,:]
    print(np.mean(magnetic_field),Vperp_squared)
    magnetic_moment = mass_spec*Vperp_squared/(2*np.mean(magnetic_field))
    plt.plot(restart_times/tcmin,magnetic_moment,color=colors[j],marker='o')
plt.show()