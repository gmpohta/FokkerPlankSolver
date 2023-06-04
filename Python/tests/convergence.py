import ctypes
import numpy as np
import matplotlib.pyplot as plt
import os
import time
import matplotlib as mpl

'''
Method convergence test, see sect 7 in Readme
'''

if os.name == 'posix':
    path = os.path.abspath(os.path.join(os.getcwd(), "../../fokker_plank_solver/fokker_plank_solver.so"))
elif os.name == 'nt':
    path = os.path.abspath(os.path.join(os.getcwd(), "../../fokker_plank_solver/fokker_plank_solver.dll"))
else:
    print('This operating system is not yet supported')
    exit()

libc=ctypes.cdll.LoadLibrary(path)

libc.get_grids.restype=None
libc.create_solver.restype=None
libc.solve.restype=None

''' Arguments for create_solver function:
    general_param = [nE, nmu, nz, nz_loop_param, switch_reverse_current, switch_turb]
    grid_param=[Emin,zmax,tmax,dt0,dt_max,gamma, r_gridE, r_gridz]
    inf_param=[n_inf, log_inf, conductivity_inf]
    loop_param=[z_loop,B_loop,n_loop,coulomb_log,conductivity]
    turb_param=[lamb_turb0, E_turb0, a_turb]
'''
libc.create_solver.argtypes=(np.ctypeslib.ndpointer(dtype=np.uint), np.ctypeslib.ndpointer(dtype=np.float64),
                             np.ctypeslib.ndpointer(dtype=np.float64), np.ctypeslib.ndpointer(dtype=np.float64),
                             np.ctypeslib.ndpointer(dtype=np.float64))

''' Arguments for solve function:
    time
    S_out =[]
    J_out =[]
    f_out=[]
    n_fast_out =[]
'''
libc.solve.argtypes=(ctypes.POINTER(ctypes.c_double),np.ctypeslib.ndpointer(dtype=np.float64), np.ctypeslib.ndpointer(dtype=np.float64),
                             np.ctypeslib.ndpointer(dtype=np.float64), np.ctypeslib.ndpointer(dtype=np.float64))

'''Arguments for get_grids function:
    E_out=[]
    mu_out=[]
    z_out=[]
    dlnB/dz=[]
    n_loop=[]
'''
libc.get_grids.argtypes=(np.ctypeslib.ndpointer(dtype=np.float64), np.ctypeslib.ndpointer(dtype=np.float64),
                        np.ctypeslib.ndpointer(dtype=np.float64),np.ctypeslib.ndpointer(dtype=np.float64), np.ctypeslib.ndpointer(dtype=np.float64))

###const###
qe=1.6e-19 #Kl
me=9.1e-31 #kg
c=3e8 #m/s
eps0=8.85e-12 #F/m

###general_param###
nE=40
nz=40
nz_loop_param=121
nmu=40
switch_turb=0
switch_reverse_current=1
general_param = np.array([nE, nmu, nz, nz_loop_param, switch_reverse_current, switch_turb],dtype=np.uint)

###grid_param###
Emin=3e-2 #MeV
zmax=6e7 # Length of the loop, m
tmax=3.0 # Time, s
r_gridz=0.25*zmax #m
r_gridE=Emin*3.0#MeV
gamma=1
dt0=tmax*1/50 # Time, s
dt_max=tmax*1/50# Time, s
grid_param=np.array([Emin,zmax,tmax,dt0,dt_max,gamma, r_gridE, r_gridz])

###inf_param###
n_inf=1e20 #m-3
log_inf=20.0
kTe_loop=100.0*qe
conductivity_inf=3.0/2.0/(2.0*np.pi)**0.5/qe/qe/me**0.5/log_inf*(4.0*np.pi*eps0)**2*kTe_loop**1.5 #Omh-1m-1
inf_param=np.array([n_inf, log_inf, conductivity_inf])

###loop_param###
z_loop=np.linspace(0,zmax,nz_loop_param) #m
B_loop=np.exp((z_loop/zmax-0.5)**2*4*np.log(5)) #G
n_loop_min=2.5e16#m-3
n_loop=n_inf-(n_inf-n_loop_min)*np.exp(-1e3*(z_loop/zmax-0.5)**16) #m-3
coulomb_log=20.0*np.ones(nz_loop_param)
conductivity=3.0/2.0/(2.0*np.pi)**0.5/qe/qe/me**0.5/coulomb_log*(4.0*np.pi*eps0)**2*kTe_loop**1.5#1/Omh/m
loop_param=np.array([z_loop,B_loop,n_loop,coulomb_log,conductivity])
loop_param=loop_param.reshape(5*nz_loop_param)

###turb_param###
a_turb=1.0
E_turb0=2e-2#MeV
lamb_turb0=1e6 #m
kminE=1.0-np.exp(-Emin/r_gridE)
E_turb=-r_gridE*np.log(1-np.linspace(kminE,1,nE+1))
turb_param=np.ones(nE+1)

###
libc.create_solver(general_param, grid_param, inf_param, loop_param, turb_param)

###get_grids###
E_grid=np.zeros(nE+1)
mu_grid=np.zeros(nmu+1)
z_grid=np.zeros(2*nz+1)
dlnB_out=np.zeros(2*nz+1)
n_loop_out=np.zeros(2*nz+1)
libc.get_grids(E_grid,mu_grid,z_grid,dlnB_out,n_loop_out)

###Initial distr function###
f=np.zeros((nE+1,nmu+1,2*nz+1))#m-3MeV-1

###Initial current###
J=np.zeros(2*nz+1) #A/m2

###output_density of fast electrons###
n_fast=np.zeros(2*nz+1)

###Injection function###
delt = 5
t0 = 0.3 #dimensionless
kt=2
E0 = Emin*qe*1e6/me/c/c  #dimensionless
I=7822.691892717618 # this is integral (16) in Readme

re_0=1/4/np.pi/eps0*qe**2/me/c**2 # m
def lamb0():
    return 1/4/np.pi/re_0**2/n_loop_out/log_inf #m

def a1(ii,jj):
    E=E_grid[ii]*qe/me/c/c*1e6
    b=(E**2+2.0*E)**0.5/(E+1.0)
    return -c*tmax/zmax*b*mu_grid[jj]

def a2(ii,jj):
    E=E_grid[ii]*qe/me/c/c*1e6
    b=(E**2+2.0*E)**0.5/(E+1.0)
    return c*tmax/zmax*(1-mu_grid[jj]**2)/2*b*dlnB_out*zmax-2*c*tmax*mu_grid[jj]/lamb0()/b**3/(E+1)**2

def aJ2(ii,jj):
    E=E_grid[ii]*qe/me/c/c*1e6
    b=(E**2+2.0*E)**0.5/(E+1.0)
    return -qe**2*n_loop_min*tmax/me/conductivity_inf*(1-mu_grid[jj]**2)/b

def a3(ii):
    E = E_grid[ii] * qe / me / c / c * 1e6
    b = (E**2 + 2.0 * E)**0.5 / (E + 1.0)
    return c*tmax/lamb0()/b

def aJ3(ii,jj):
    E = E_grid[ii] * qe / me / c / c * 1e6
    b = (E**2 + 2.0 * E)**0.5 / (E + 1.0)
    return -qe**2*n_loop_min*tmax/me/conductivity_inf*b*mu_grid[jj]

def a4(ii,jj):
    E = E_grid[ii] * qe / me / c / c * 1e6
    b = (E**2 + 2.0 * E)**0.5 / (E + 1.0)
    return c*tmax/lamb0()*(1-mu_grid[jj]**2)/b**3/(E+1)**2

def a5(ii,jj):
    E = E_grid[ii] * qe / me / c / c * 1e6
    b = (E**2 + 2.0 * E)**0.5 / (E + 1.0)
    return -c*tmax/lamb0()/(E**2+2*E)**1.5

def source(time):
    S=np.zeros((nE+1,nmu+1,2*nz+1)) #m-3MeV-1s-1
    t_d=time / tmax
    ft =0
    dft=0
    if t_d < t0:
        ft = t_d
        dft = 1
    else:
        ft = t0 * np.exp(-kt*(t_d - t0)**2)
        dft = -kt * 2*(t_d  - t0) * ft

    fz = np.zeros(2*nz+1)
    dfz = np.zeros(2*nz+1)

    for kk in  range(2*nz+1):
        z_d=z_grid[kk]/zmax
        if z_d>=0 and z_d<=1:
            fz[kk] = np.sin(np.pi * z_d)
            dfz[kk]=np.pi*np.cos(np.pi*z_d)
        else:
            fz[kk]=0
            dfz[kk]=0

    for ii in range(nE + 1):
        E_d=E_grid[ii] * qe * 1e6 / me / c / c
        fE=E0**delt/E_d**delt
        dfE = -delt * E0**delt / E_d**(delt + 1)
        for jj in range(nmu + 1):
            fmu = 1/10  - mu_grid[jj]**3/10
            dfmu = -3 * mu_grid[jj]**2/10
            d2fmu = -6 * mu_grid[jj]/10

            Ja = 2*fz * ft / 5.0/10 * E0**delt * I
            S[ii,jj,:]=dft*fE*fmu*fz-a1(ii,jj)*ft*fE*fmu*dfz-a5(ii,jj)*ft*fE*fmu*fz-\
                (a2(ii,jj)+aJ2(ii,jj)*Ja)*ft*fE*dfmu*fz-(a3(ii)+aJ3(ii,jj)*Ja)*ft*dfE*fmu*fz-a4(ii,jj)*ft*fE*d2fmu*fz
            S[ii, jj, :] *=n_loop_min*qe*1e6/me/c/c/tmax #convert to MeV-1m-3s-1
    return S

###calculation results for each point in time###
f_out=[]
J_out=[]
n_fast_out=[]
t=[]
tmp_t = ctypes.c_double(0.0)
f_out.append(f.copy())
J_out.append(J.copy())
t.append(tmp_t.value)
ct=1

###solution cycle##
time._startTime = time.time()
while tmp_t.value<tmax:
    S=source(tmp_t.value)
    libc.solve(tmp_t,S,J,f,n_fast)
    t.append(tmp_t.value)
    f_out.append(f.copy())
    J_out.append(J.copy())
    n_fast_out.append(n_fast.copy())
    print(tmp_t.value,ct)
    ct+=1
print('Elapsed time: ', time.time() - time._startTime)

f_out=np.array(f_out)
J_out=np.array(J_out)
n_fast_out=np.array(n_fast_out)
t=np.array(t)

def analitic_solve(t,E,m,z):
    nt=len(t)
    fa=np.zeros((nt,nE+1,nmu+1,2*nz+1))
    Ja=np.zeros((nt,2*nz+1))
    fz=np.zeros(2*nz+1)
    for kk in range(2 * nz + 1):
        z_d=z[kk] / zmax
        if z_d >= 0 and z_d <= 1:
            fz[kk] = np.sin(np.pi * z_d)
        else:
            fz[kk] = 0

    for ct in range(nt):
        t_d=t[ct] / tmax
        if  t_d>= 0 and t_d <= t0:
            ft = t_d
        else:
            ft = t0 * np.exp(-kt*(t_d - t0)**2)
        for jj in range(nmu+1):
            fm=1/10-m[jj]**3/10
            for ii in range(nE + 1):
                E_d=E[ii]*qe*1e6/me/c/c
                fE=(E0/E_d) ** delt
                fa[ct,ii,jj,:]=fE*fm*fz*ft*n_loop_min*qe*1e6/me/c/c #convert to MeV-1m-3

        Ja[ct,:]=2*fz*ft/5/10*E0**delt*I*qe*c*n_loop_min #convert to A/m2
    return fa ,Ja
##
fa,Ja=analitic_solve(t,E_grid,mu_grid,z_grid)

#find indices which correspond max error
ind=np.unravel_index(np.argmax(abs(fa-f_out)),fa.shape)
mi=ind[2]
zi=ind[3]
Ei=ind[1]
ti=ind[0]
print('error',abs(fa[ind]-f_out[ind])/(np.max(abs(fa))))
print(Ei,E_grid[Ei],mi,mu_grid[mi],zi,z_grid[zi])

front_size_tick=16
front_size_labels=16
def plot_3D(x,y,z1,z2,labels,legend,zscale='normal'):
    fig=plt.figure()
    axes=fig.add_subplot(projection='3d')
    axes.plot_wireframe(x,y,z1,color='black',label=legend[0])
    axes.plot_wireframe(x, y, z2, color='red',label=legend[1])
    axes.legend()
    plt.xticks(fontsize=front_size_tick)
    plt.yticks(fontsize=front_size_tick)
    for t_ax in axes.zaxis.get_major_ticks(): t_ax.label.set_fontsize(front_size_tick)
    axes.set_xlabel(labels['x'],fontsize=front_size_labels,labelpad=10)
    axes.set_ylabel(labels['y'],fontsize=front_size_labels,labelpad=10)
    axes.zaxis.set_rotate_label(False)
    axes.set_zlabel(labels['z'],fontsize=front_size_labels,rotation=90,labelpad=15)

    def log_tick_formatter(val, pos=None):
        return "$10^{{{0:.0f}}}$".format(val)

    if zscale=='log':
        axes.zaxis.set_major_formatter(mpl.ticker.FuncFormatter(log_tick_formatter))
        axes.zaxis.set_major_locator(mpl.ticker.MaxNLocator(integer=True))
    return axes

##
tg,Eg=np.meshgrid(E_grid,t)
labels={'x':'Energy, MeV','y':'Time,s','z':'$\mathit{f(E,\mu,z,t)}, MeV^{-1}m^{-3}$'}
leg=['numerical','analytical']
ax=plot_3D(tg,Eg,np.log10(f_out[:,:,mi,zi]),np.log10(fa[:,:,mi,zi]),labels,leg,zscale='log')

##
tg,mg=np.meshgrid(mu_grid,t)
labels={'x':'Cos pith angle','y':'Time,s','z':'$\mathit{f(E,\mu,z,t)}, MeV^{-1}m^{-3}$'}
ax=plot_3D(tg,mg,f_out[:,Ei,:,zi],fa[:,Ei,:,zi],labels,leg)

##
tg,zg=np.meshgrid(z_grid,t)
labels={'x':'Distance, m','y':'Time,s','z':'$\mathit{f(E,\mu,z,t)}, MeV^{-1}m^{-3}$'}
ax=plot_3D(tg,zg,f_out[:,Ei,mi,:],fa[:,Ei,mi,:],labels,leg)

##
labels={'x':'Distance, m','y':'Time,s','z':'$\mathit{J(z,t)}, A/m^2$'}
ax=plot_3D(tg,zg,J_out,Ja,labels,leg)

plt.show()
