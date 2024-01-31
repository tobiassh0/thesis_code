
def jld_to_pkl(loc='',frac=1):
    f = h5py.File(loc+'/solutions2D_.jld',"r")
    keys=f.keys()
    # for k in keys: print(k)
    w0 = f["w0"][()]
    k0 = f["k0"][()]
    sols = f["plasmasols"][()]
    sshape = sols.shape[0]
    w = np.zeros(int(sshape/frac+1)) ; dw = np.zeros(int(sshape/frac+1))
    kpara = np.zeros(int(sshape/frac+1)) ; kperp = np.zeros(int(sshape/frac+1))
    i=0
    for item in sols[::frac]:
        print(i)
        w[i],dw[i] = f[item][()][0]
        kpara[i],kperp[i] = f[item][()][1]
        i+=1
    dumpfiles(w,'frequency')
    dumpfiles(dw,'growthrates')
    dumpfiles(kpara,'parallelwavenumber')
    dumpfiles(kperp,'perpendicularwavenumber')
    dumpfiles([w0,k0],'w0k0')
    return None

def read_all_data(loc=''):
    home = os.getcwd()
    os.chdir(loc) # change to where data is saved
    # load all relevant data
    w0,k0=read_pkl('w0k0')
    w=read_pkl('frequency')
    dw=read_pkl('growthrates')
    kpara=read_pkl('parallelwavenumber')
    kperp=read_pkl('perpendicularwavenumber')
    os.chdir(home) # change back to home
    return w0,k0,w,dw,kpara,kperp

def plot_growth_k2d(kpara,kperp,dw,norm=[None,None],cmap='viridis'):
    # cm = plt.cm.get_cmap('viridis')
    # sc = plt.scatter(kperp/w0,kpara/k0,c=dw/w0,cmap=cm,s=1)
    # plt.colorbar(sc)
    # plt.show()
    x=np.unique(kperp)
    y=np.unique(kpara)
    X,Y=np.meshgrid(x,y)
    Z=dw.reshape(len(y),len(x))
    plt.colormesh(X,Y,Z/w0,cmap=cm)
    plt.show()
    return None

def plot_freq_kpara(w,dw,kpara,norm=[None,None],colors='viridis'):
    return None

def plot_freq_kperp(w,dw,kperp,norm=[None,None],colors='viridis'):
    return None

def plot_freq_growth(w,dw,kpara,norm=[None,None],colors='viridis'):
    thresh = dw > 0
    plt.scatter(w[thresh]/w0,dw[thresh]/w0,c=cmap)
    plt.savefig('freq_growth.pdf')
    return None

def make_all_plots(alldata=None,cmap='viridis'):
    w0,k0,w,dw,kpara,kperp = alldata
    plot_growth_k2d(kpara,kperp,dw,norm=[w0,k0],cmap=cmap)
    # plot_freq_kpara(w,dw,kpara,norm=[w0,k0],colors=cmap)
    # plot_freq_kperp(w,dw,kperp,norm=[w0,k0],colors=cmap)
    # plot_freq_growth(w,dw,kpara,norm=[w0,k0],colors=cmap)
    return None


if __name__ == '__main__':
    from func_load import *
    import h5py
    solloc = "/home/space/phrmsf/Documents/ICE_DS/JET26148/default_params_with_Triton_concentration/run_2.07_0.11_-0.646_0.01_0.01_15.0_3.5__1.0_4.0_1.7e19_0.00015_1024"
    jld_to_pkl(loc=solloc)
    data=read_all_data()
    make_all_plots(alldata=data)
