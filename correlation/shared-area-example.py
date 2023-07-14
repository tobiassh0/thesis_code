from list_new import *
from scipy import interpolate

# shared area example
t = np.linspace(0,10,1000)
f = 2*const.PI

peak_area = []
peak_spline=[]
N=1e4
vals = [0.,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5,2.75,3.0,3.25,3.5]
peak_area=np.zeros((len(vals),int(N)))
peak_spline=np.zeros((len(vals),int(N)))
## Monte Carlo method
for j in range(len(vals)):
	## signals
	a1 = 1. ; b1 = 1
	a2 = 6. ; b2 = b1
	noise_floor = vals[j]*b1
	pa = np.zeros(int(N)) ; ps = np.zeros(int(N))
	for i in range(int(N)):
		noise1 = np.random.normal(0,noise_floor,len(t)) + 5
		noise2 = np.random.normal(0,noise_floor,len(t)) + 5
		sig1 = np.exp(-b1*(t-a1)**2)+noise1
		sig2 = np.exp(-b2*(t-a2)**2)+noise2

		## shared area between signals
		areas = shared_area(sig1,sig2) # from list_new

		## fitting spline
		knot_numbers = 20
		t_new = np.linspace(0, 1, knot_numbers+2)[1:-1]
		q_knots = np.quantile(t, t_new)
		tf,c,k = interpolate.splrep(t, areas, t=q_knots, s=1)
		yfit = interpolate.BSpline(tf,c,k)(t)
		## summing total
		pa[i]=t[np.argmax(areas)]
		ps[i]=t[np.argmax(yfit)]
	peak_area[j,] = (pa)
	peak_spline[j,] = (ps)

mean_area=np.mean(peak_area,axis=1)
mean_spline=np.mean(peak_spline,axis=1)
std_area=np.std(peak_area,axis=1)
std_spline=np.std(peak_spline,axis=1)
plt.errorbar(vals,mean_area,yerr=std_area,fmt='o',color='b')
plt.errorbar(vals,mean_spline,yerr=std_spline,fmt='s',color='g')
plt.plot(vals,mean_area,'-o',vals,mean_spline,'--s')

plt.gca().ticklabel_format(useOffset=False)
plt.ylabel('Offset, '+r'$\tau$',**tnrfont)
plt.xlabel('Noise amp./Signal amp.',**tnrfont)
plt.savefig('noisy-shared-area_MC_errorbars.png')
plt.show() ; sys.exit()

