
from func_load import *

home = '/storage/space2/phrmsf/lowres_D_He3/'
os.chdir(home)
sims = np.sort([i for i in os.listdir() if 'p_90' in i])
xiHe3 = [int(i[2:4]) for i in sims]
energy_max=[]
energy_integral=[]
minspec = 'Protons'
for sim in sims:
    simloc = getSimulation(sim)
    times = read_pkl('times')
    dt = (times[-1]-times[0])/len(times)
    energy_minspec = read_pkl(minspec+'_KEdens')
    energy_minspec = energy_minspec - np.mean(energy_minspec[0:10]) # change in energy dens
    energy_abs = np.abs(energy_minspec)
    # plt.plot(times,energy_abs)
    energy_max.append(np.max(energy_abs))
    energy_integral.append(np.sum(energy_minspec*dt))
    os.chdir(home)
# plt.show()
plt.scatter(xiHe3[0],1)
for i in range(1,len(energy_integral)):
    plt.scatter(xiHe3[i]/100,(energy_integral[i]-energy_integral[0])/energy_integral[0])
plt.show()