import ctypes
import numpy as np
import os
import matplotlib.pyplot as plt
import time

'''
Turbulent diffusion test, see sect 8 in Readme
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
dt_max=tmax*0.02# Time, s
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
        for kk in range(2*nz+1):
            f[ii,jj,kk]=10.0*np.exp(-5e2*(mu_grid[jj]-0.7)**2)*(Emin/E_grid[ii])**delta*np.exp(-(z_grid[kk]/zmax-0.5)**2/s0/s0)*n_loop_min*qe*1e6/me/c/c #m-3MeV-1

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


def plot_3D(x,y,z,title,labels):
    front_size_tick=16
    front_size_labels=16
    fig=plt.figure()

    axes=fig.add_subplot(projection='3d')
    axes.set_title(title)

    axes.plot_wireframe(x,y,z,color='black')

    plt.xticks(fontsize=front_size_tick)
    plt.yticks(fontsize=front_size_tick)
    axes.zaxis.set_tick_params(labelsize=front_size_tick)

    axes.set_xlabel(labels['x'],fontsize=front_size_labels,labelpad=10)
    axes.set_ylabel(labels['y'],fontsize=front_size_labels,labelpad=10)
    axes.zaxis.set_rotate_label(False)
    axes.set_zlabel(labels['z'],fontsize=front_size_labels,rotation=90,labelpad=15)
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
    #axes.clabel(cs,cs.levels[::2],colors="black",fontsize=16, inline=1)
    manual_locations = [(0.3, 1e7), (0.8, 1.5e7), (0.28, 1.95e7), (0.25, 2.4e7), (0.32, 2.7e7), (0.3, 5.1e7),(0.8,4.5e7),(0.28,4.2e7)]
    axes.clabel(cs, colors="black",fontsize=16, manual=manual_locations)
    axes.tick_params(axis='both',which='major',labelsize=front_size_tick)
    axes.grid()
    return axes

plt.plot(z_grid,dlnB_out)
plt.plot(z_grid,dlnB_out,'or')
plt.xlabel('Distance along loop, m')
plt.ylabel('$\dfrac{d\ln B}{dz}, m^{-1}$')
plt.grid()
plt.figure()
plt.plot(z_grid,n_loop_out)
plt.plot(z_grid,n_loop_out,'or')
plt.xlabel('Distance along loop, m')
plt.ylabel('Density thermal plasma, $m^{-3}$')
plt.grid()

tg,zg=np.meshgrid(z_grid,t[1:])
labels={'z':'$\mathit{n}, m^{-3}$','x':'Time,s ','y':'Distance, m'}
ax=plot_3D(zg,tg,n_fast_out,'Distance',labels)

title='$\mathit{n, rel.units}$'
labels={'y':'Distance, m','x':'Time, s'}
lev=np.concatenate((np.linspace(1e-4,1e-2,5),np.linspace(1.5e-2,1,20)))
ax=plot_contour(zg[:,1:-1],tg[:,1:-1],n_fast_out[:,1:-1]/np.max(n_fast_out[:,1:-1]),title,labels,lev)
ax.set_ylim([-zmax/20,zmax+zmax/20])
plt.show()
