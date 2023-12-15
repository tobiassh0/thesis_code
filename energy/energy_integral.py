
from func_load import *

home = '/storage/space2/phrmsf/lowres_D_He3/'
os.chdir(home)
sims = np.sort([i for i in os.listdir() if 'p_90' in i])[1:] # ignoring 0% case
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
# scatter plot each point w.r.t the 5% case
plt.scatter(xiHe3[0]/100,0)
for i in range(1,len(energy_integral)):
    plt.scatter(xiHe3[i]/100,(energy_integral[i]-energy_integral[0])/energy_integral[0])
# plot straight 1:1 line
plt.plot([0,0.5],[0,0.5],linestyle='--',color='darkgrey')
# calculate ODR fit with errors
int_diff = [0] + [(energy_integral[i] - energy_integral[0])/energy_integral[0] for i in range(1,len(energy_integral))]
xiHe3 = np.array(xiHe3)
xiarr = np.arange(0,0.5,0.01)
params,params_err= ODR_fit(xiHe3/100,int_diff,sx=None,sy=None,beta0=[1,0],curve='linear')
plt.plot(xiarr,xiarr*params[0]+params[1],color='r',linestyle='--')
print(params,params_err) # [ 0.98460096 -0.01862662] [0.16437624 0.04512366]
plt.ylabel(r'$(S_{\xi_{He3}}-S_{5\%})/S_{5\%}$',**tnrfont)
plt.xlabel(r'$\xi_{He3}$',**tnrfont)
plt.ylim(-0.05,0.5) ; plt.xlim(0,0.5)
plt.show()

"""
Find that this is basically a straight line, therefore the change in the total change in energy density of
the energetic protons is directly proportional to the quantity of He3 present in the plasma.
"""