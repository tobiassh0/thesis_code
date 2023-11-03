import os

#rootdir = '/storage/space2/phrmsf/lowres_D_He3/'
rootdir = '/home/space/phrmsf/Documents/thesis_code/'

filename = '_KE'
newfilename = '_KEdens'

for subdir, dirs, files in os.walk(rootdir):
	for f in files:
		if filename in f:
			print(os.path.join(subdir,f))
			os.rename(os.path.join(subdir,f),os.path.join(subdir,str(f).replace(filename,newfilename)))
			print(f)
