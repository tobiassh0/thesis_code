from func_load import *

os.chdir('/storage/space2/phrmsf/lowres_D_He3/')
sims = [i for i in os.listdir() if 'p_90' in i]
hlabels = np.sort([int(i[2:4]) for i in sims])
xlabel=r'$[\Delta u_{He3}/\Delta u_D]_{max}$'
ylabel=r'$(\xi_{He3}/\xi_{D})(m_{D}/m_{He3})(q_{He3}/q_D)^2$'
majIons_edens_ratio(sims,species=['He3','Deuterons'],time_norm=r'$\tau_{cD}$',labels=hlabels,xlabel=xlabel,ylabel=ylabel,lims=((0,5),(0,5)))

