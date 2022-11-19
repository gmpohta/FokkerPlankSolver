import ctypes
import numpy as np
import matplotlib.pyplot as plt
import os
import time
from scipy.integrate import quad

'''
Compare with reduced equation, see sect 9 in Readme
'''

libc=ctypes.cdll.LoadLibrary(os.getcwd()+'//FokkerPlankSolver.dll')

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
nz=50
nz_loop_param=121
nmu=40
switch_turb=1
switch_reverse_current=1
general_param = np.array([nE, nmu, nz, nz_loop_param, switch_reverse_current, switch_turb],dtype=np.uint)

###grid_param###
Emin=3e-2 #MeV
zmax=6e7 # Length of the loop, m
tmax=1.0 # Time, s
r_gridz=0.25*zmax #m
r_gridE=Emin*5.0#MeV
gamma=0.1
dt0=tmax*1e-4 # Time, s
dt_max=tmax*0.005# Time, s
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
n_loop=n_inf-(n_inf-n_loop_min)*np.exp(-5e3*(z_loop/zmax-0.5)**16) #m-3
coulomb_log=20.0*np.ones(nz_loop_param)
conductivity=3.0/2.0/(2.0*np.pi)**0.5/qe/qe/me**0.5/coulomb_log*(4.0*np.pi*eps0)**2*kTe_loop**1.5#1/Omh/m
loop_param=np.array([z_loop,B_loop,n_loop,coulomb_log,conductivity])
loop_param=loop_param.reshape(5*nz_loop_param)

###turb_param###
a_turb=1.0
lamb_turb0=1e6 #m
kminE=1.0-np.exp(-Emin/r_gridE);
E_turb=-r_gridE*np.log(1-np.linspace(kminE,1,nE+1))
turb_param=lamb_turb0*(Emin/E_turb)**a_turb

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
s0=3e6/zmax
delta=5.0
f=np.zeros((nE+1,nmu+1,2*nz+1))
for ii in range(nE):
    for jj in range(nmu):
        f[ii,jj,:]=10.0*np.exp(-5e2*(mu_grid[jj]-0.7)**2)*(Emin/E_grid[ii])**delta*np.exp(-(z_grid/zmax-0.5)**2/s0/s0)*n_loop_min*qe*1e6/me/c/c #m-3MeV-1

###Initial current###
J=np.zeros(2*nz+1) #A/m2

###output_density of fast electrons###
n_fast=np.zeros(2*nz+1)

###Injection function###
S=np.zeros((nE+1,nmu+1,2*nz+1))

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

def simple_equation(n_loop,turb_param,Emin_in,rz_in,rE_in,nz,nE):
    r0_e = 1.0/4.0/np.pi/eps0*qe*qe/me/c/c
    lambda0_e = 1.0/4.0/np.pi/r0_e/r0_e

    # number od intervals
    nt = 200

    # len of intervals
    dz = 2/nz
    dt = 1/nt
    # arrays of solution
    f = np.nan*np.ones((nE + 1, 2 * nz + 1, nt + 1))
    I = np.nan*np.ones((nE + 1, 2 * nz + 1, nt + 1))

    # grid t (dimensionless)
    t = np.linspace(0, 1, nt + 1)

    # Non uniform grids
    # for E
    Emin=Emin_in/me/c/c*qe*1e6 #dimensionless
    rE = rE_in*qe/me/c/c*1e6 #dimensionless
    kmin = 1 - np.exp(-Emin / rE)
    dkE = (1 - kmin) / nE
    kE = np.linspace(kmin, 1, nE + 1)
    E = -rE * np.log(1 - kE[:-1])
    dkdE = rE / (1 - kE)
    # for z
    rz = rz_in/zmax #dimensionless
    dkz = 1 / nz
    kz2 = np.linspace(0, 1, nz + 1)
    kz1 = np.linspace(-1, 0, nz + 1)
    z = np.concatenate((rz * np.log(1 + kz1[:-1]), -rz * np.log(1 - kz2)))
    dZdK = np.concatenate((rz / (1 + kz1[:-1]), rz / (1 - kz2)))

    # Initial
    def f_mu(mu):
        return np.exp(-5e2 * (mu - 0.7) ** 2)
    int_f_mu= quad(f_mu,-1, 1)[0]
    for ii in range(nE):
        f[ii, :, 0] = 10*(Emin / E[ii]) ** delta * np.exp(-z ** 2 / s0 ** 2) *int_f_mu  # initial f, dimensionless
    f[:, -1, :] = 0  # z
    f[:, 0, :] = 0  # z
    f[-1, :, :] = 0  # E

    def f_lam0(kk):
        lam0 = lambda0_e/ n_loop[kk] / log_inf # m
        return lam0

    def kE4(ii, kk):
        b = (E[ii] * (E[ii] + 2)) ** 0.5 / (E[ii] + 1)
        lam0 = f_lam0(kk)
        return c * tmax / lam0 / b

    def kE2(ii, kk):
        lam0 = f_lam0(kk)
        return -c * tmax / lam0 / (E[ii] ** 2 + 2 * E[ii]) ** 1.5

    def kturb(ii):
        b = (E[ii] * (E[ii] + 2)) ** 0.5 / (E[ii] + 1)
        return c * tmax / zmax ** 2 * turb_param[ii]* b / 3

    def tridiag_alg(a, b, c, r):
        # b - main diagonal
        # a - lower diaganal
        # c - upper diaganal
        n = len(b)
        alf = np.zeros(n)
        beta = np.zeros(n)
        out = np.zeros(n)
        alf[0] = -c[0] / b[0]
        beta[0] = r[0] / b[0]
        for ii in range(1, n):
            alf[ii] = -c[ii] / (b[ii] + alf[ii - 1] * a[ii])
            beta[ii] = (r[ii] - a[ii] * beta[ii - 1]) / (b[ii] + alf[ii - 1] * a[ii])
        out[-1] = beta[-1]
        for ii in reversed(range(n - 1)):
            out[ii] = alf[ii] * out[ii + 1] + beta[ii]
        return out

    for ct in range(nt):
        f1 = np.ones((nE + 1, 2 * nz + 1)) * np.nan
        f1[-1, :] = 0  # E
        f1[:, -1] = 0  # z
        f1[:, 0] = 0  # z
        for ii in reversed(range(nE)):
            a4 = kE4(ii, range(2*nz+1))
            a2 = kE2(ii, range(2*nz+1))
            aturb = kturb(ii)

            op_E = a4 / dkE / dkdE[ii] * (f[ii + 1, :, ct] - f[ii, :, ct])

            op_z = np.ones(2 * nz + 1) * np.nan
            for kk in range(1, 2 * nz):
                op_z[kk] = aturb / dkz ** 2 / dZdK[kk] ** 2 * (
                            f[ii, kk + 1, ct] - 2 * f[ii, kk, ct] + f[ii, kk - 1, ct])
            sou = op_E + op_z + a2 * f[ii, :, ct]  ########
            f1[ii, :] = (sou + a4 * dt / dkE / dkdE[ii] * f1[ii + 1, :]) / (1 + a4 * dt / dkE / dkdE[ii] - a2 * dt)

        for ii in range(nE):
            aturb = kturb(ii)
            ay = -aturb * dt / dkz ** 2 / dZdK ** 2 / 2
            cy = -aturb * dt / dkz ** 2 / dZdK ** 2 / 2
            ry = np.ones(2 * nz + 1) * np.nan
            for kk in range(1, 2 * nz):
                ry[kk] = f1[ii, kk] * dt - aturb * dt / 2 / dkz ** 2 / dZdK[kk] ** 2 * (
                            f[ii, kk + 1, ct] - 2 * f[ii, kk, ct] + f[ii, kk - 1, ct]) + f[ii, kk, ct]
            by = 1 + aturb * dt / dkz ** 2 / dZdK ** 2
            f[ii, 1:-1, ct + 1] = tridiag_alg(ay[1:-1], by[1:-1], cy[1:-1], ry[1:-1])
        print(ct, ' in ', nt - 1)

    n_fast = 0
    for ii in range(1, nE):
        n_fast += (f[ii, :, :] + f[ii - 1, :, :]) * dkdE[ii] * (kE[ii] - kE[ii - 1]) / 2

    z = (z+1/2) * zmax #m
    t *= tmax #s
    n_fast *= n_loop_min #m-3

    return n_fast, z, t

time._startTime = time.time()
n_fast_simpl,z_simpl,t_simpl=simple_equation(n_loop_out,turb_param,Emin,r_gridz,r_gridE,nz,nE)
print('Elapsed time: ', time.time() - time._startTime)

def plot_3D(x,y,z,x2,y2,z2,title,labels,legend):
    front_size_tick=16
    front_size_labels=16
    fig=plt.figure()
    
    axes=fig.add_subplot(projection='3d')
    axes.set_title(title)

    axes.plot_wireframe(x,y,z,color='black',label=legend[0])
    axes.plot_wireframe(x2, y2, z2, color='red', label=legend[1])

    plt.xticks(fontsize=front_size_tick)
    plt.yticks(fontsize=front_size_tick)
    for t_ax in axes.zaxis.get_major_ticks(): t_ax.label.set_fontsize(front_size_tick)
    axes.set_xlabel(labels['x'],fontsize=front_size_labels,labelpad=10)
    axes.set_ylabel(labels['y'],fontsize=front_size_labels,labelpad=10)
    axes.zaxis.set_rotate_label(False)
    axes.set_zlabel(labels['z'],fontsize=front_size_labels,rotation=90,labelpad=15)
    axes.legend()
    def mjrFormatter(x,pos):
        return "$10^{{{0:.0f}}}$".format(x)
    return axes

def plot_contour(x,y,z,title,labels,levels):
    front_size_tick=20
    front_size_labels=20

    fig=plt.figure()
    axes=plt.axes()
    axes.set_xlabel(labels['x'],fontsize=front_size_labels)
    axes.set_ylabel(labels['y'],fontsize=front_size_labels)
    axes.set_title(title,fontsize=front_size_labels)
    cs=axes.contour(x,y,z,levels)
    axes.clabel(cs,cs.levels[::2],colors="black",fontsize=16)
    axes.tick_params(axis='both',which='major',labelsize=front_size_tick)
    axes.grid()
    return axes

tg,zg=np.meshgrid(z_grid,t[1:])
tg_simpl,zg_simpl=np.meshgrid(z_simpl,t_simpl)
n_fast_simpl=n_fast_simpl.swapaxes(0, 1)

labels={'z':'$\mathit{n(s,t)}, m^{-3}$','x':'Time,s','y':'Distance, m'}
leg=['full equation','reduced equation']
axes=plot_3D(zg,tg,n_fast_out,zg_simpl,tg_simpl,n_fast_simpl,'n_fast',labels,leg)

title='$\mathit{n, rel.units}$'
labels={'y':'Distance, m','x':'Time, s'}
lev=np.concatenate((np.linspace(1e-4,1e-2,5),np.linspace(1.5e-2,1,20)))
ax=plot_contour(zg_simpl[:,1:-1],tg_simpl[:,1:-1],n_fast_simpl[:,1:-1]/np.max(n_fast_simpl[:,1:-1]),title,labels,lev)
ax.set_ylim([-zmax/20,zmax+zmax/20])

title='$\mathit{n, rel.units}$'
labels={'y':'Distance, m','x':'Time, s'}
lev=np.concatenate((np.linspace(1e-4,1e-2,5),np.linspace(1.5e-2,1,20)))
ax=plot_contour(zg[:,1:-1],tg[:,1:-1],n_fast_out[:,1:-1]/np.max(n_fast_out[:,1:-1]),title,labels,lev)
ax.set_ylim([-zmax/20,zmax+zmax/20])

plt.show()