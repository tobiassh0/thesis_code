from func_load import *

os.chdir('/storage/space2/phrmsf/lowres_D_He3/')
sims = [i for i in os.listdir() if 'p_90' in i]
majIons_edens_ratio(sims,species=['He3','Deuterons'],time_norm=r'$\tau_{cD}$',labels=[15,22,25])

