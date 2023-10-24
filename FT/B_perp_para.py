
from func_load import *

simloc = getSimulation('/storage/space2/phrmsf/ECRH/ECRH_JT60U_5_2')
Bz = load_batch_fieldmatrix([],'Magnetic_Field_Bz')
Bx = load_batch_fieldmatrix([],'Magnetic_Field_Bx')
By = load_batch_fieldmatrix([],'Magnetic_Field_By')
# try loading para
Bpara = read_pkl('Magnetic_Field_Bpara')
Bperp = read_pkl('Magnetic_Field_Bperp')

## change in magnetic fields
dBz = Bz - np.mean(Bz[0:10,:])
dBx = Bx - np.mean(Bx[0:10,:])
dBy = By - np.mean(By[0:10,:])

theta = 86.3*const.PI/180. # deg
# rotation is +ve for counter-clockwise in a R-handed cartesian coordinate system
phi = -(const.PI/2 - theta)
cos = np.cos(phi) ; sin = np.sin(phi)
Bperp = np.sqrt((cos**2)*(Bx**2) + (sin**2)*(Bz**2) + (By**2)) 
Bpara = np.sqrt((sin**2)*(Bx**2) + (cos**2)*(Bz**2))
dumpfiles(Bperp,'Magnetic_Field_Bperp')
dumpfiles(Bpara,'Magnetic_Field_Bpara')

dBpara = Bpara - np.mean(Bpara[0:10,:])
dBperp = Bperp - np.mean(Bperp[0:10,:])

clim = (-4,6)
fig,axs = plt.subplots(nrows=3,ncols=2,figsize=(10,10))
## z
im00 = axs[0,0].imshow(np.log10(dBz),**kwargs,clim=clim)
axs[0,0].set_title('dBz')
## para
im01 = axs[0,1].imshow(np.log10(dBpara),**kwargs,clim=clim)
axs[0,1].set_title('dBpara')
## x
im10 = axs[1,0].imshow(np.log10(dBx),**kwargs,clim=clim)
axs[1,0].set_title('dBx')
## perp
im11 = axs[1,1].imshow(np.log10(dBperp),**kwargs,clim=clim)
axs[1,1].set_title('dBperp')
## y
im20 = axs[2,0].imshow(np.log10(dBy),**kwargs,clim=clim)
axs[2,0].set_title('dBy')

## save plot
os.chdir('/home/space/phrmsf/Documents/thesis_code/')
plt.savefig('Brot.png')
#plt.show()
plt.clf()

## FT2d
FTpara = get2dTransform(dBpara) 
FTperp = get2dTransform(dBperp)
dumpfiles(FTpara,'FT_2d_Magnetic_Field_Bpara')
dumpfiles(FTperp,'FT_2d_Magnetic_Field_Bperp')
fig,axs = plt.subplots(nrows=1,ncols=2,figsize=(10,10))
im00 = axs[0].imshow(np.log10(FTpara),**kwargs)
axs[0].set_title('FT2d dBpara')
im00 = axs[1].imshow(np.log10(FTperp),**kwargs)
axs[1].set_title('FT2d dBperp')

os.chdir('/home/space/phrmsf/Documents/thesis_code/')
plt.savefig('FT2d_Brot.png')
#plt.show()
