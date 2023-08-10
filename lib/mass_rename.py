import os

rootdir = '/home/space/phrmsf/Documents/thesis_code/'
#rootdir = '/storage/space2/phrmsf/'

for subdir, dirs, files in os.walk(rootdir):
	for f in files:
		if '_KE' in f:
			print(os.path.join(subdir,f))
			os.rename(os.path.join(subdir,f),os.path.join(subdir,str(f).replace('_KE','_KEdens')))
			print(f)
