
from func_load import *

os.chdir('/storage/space2/phrmsf/lowres_D_He3/')
sims = [i for i in os.listdir() if 'p_90' in i]
hlabels = np.sort([int(i[2:4]) for i in sims])
ylabel=r'$[\Delta u_{D}/\Delta u_{He3}]_{max}$'
xlabel=r'$(\xi_{D}/\xi_{He3})(m_{He3}/m_{D})(q_{D}/q_{He3})^2$'
majIons_edens_ratio(sims,species=['Deuterons','He3'],time_norm=r'$\tau_{cD}$',labels=hlabels,xlabel=xlabel,ylabel=ylabel,lims=((0,10),(0,10)))

