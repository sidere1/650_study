import sys
sys.path.append('.')
from opf650 import *

import os

import numpy as np
import time

# physical data

u_wind = np.array([7.0710 , 7.0710]) # 45.0 	 max upwind
# u_wind = np.array([5.5557 , 8.3146]) # 56.25 	 
# u_wind = np.array([3.8268 , 9.2387]) # 67.5 	 
# u_wind = np.array([1.9509 , 9.8078]) # 78.75 	 
# u_wind = np.array([0      , 10.0  ]) # 90.0 	 
# u_wind = np.array([-1.950 , 9.8078]) # 101.25 	 
# u_wind = np.array([-3.826 , 9.2387]) # 112.5 	 
# u_wind = np.array([-5.555 , 8.3146]) # 123.75 	 
# u_wind = np.array([-7.071 , 7.0710]) # 135.0 	 
# u_wind = np.array([-10 , 0]) # 180    pure downwind

boat_length=6.5

# parameters

drag_file = '650_m1/postProcessing/forces/0/force.dat'
controlDict_0_file = '650_basecase/system/controlDict.m0'
controlDict_1_file = '650_basecase/system/controlDict.m1'
u_file = '650_basecase/my0/U'
opf_script = './go.sh'
max_iter = 10
threshold = 0.01 # convergence criterion
relax_factor_min = 0.2
relax_factor_max = 0.6
M0_ftt_length = 15 # initial run duration, to get an initial condition
M1_ftt_length = 75 # run full duration
M1_discarded_ftt_length = 50 # run discarded initialisation duration
M0_ftt_length = 5 # initial run duration, to get an initial condition
M1_ftt_length = 25 # run full duration
M1_discarded_ftt_length = 15 # run discarded initialisation duration

# processing

for filename in os.listdir("logs"):
    file_path = os.path.join("logs", filename)
    try:
        if os.path.isfile(file_path):
            os.remove(file_path)
    except Exception as e:
        print(f'Failed to delete {file_path}. Reason: {e}')


converged = False
u_boat = np.zeros(max_iter+1)
T_all = np.zeros(max_iter+1)
D_all = np.zeros(max_iter+1)
# u_boat[0] = 7.4
# u_boat[0] = 6.8
u_boat[0] = estimate_initial_velocity(u_wind, verbose=False, plot_res=False)
relax_factors = np.linspace(0,1,max_iter)*(relax_factor_max-relax_factor_min)+relax_factor_min

t0 = time.time()
t1 = time.time()
for i_iter in np.arange(max_iter):
    print(f'Iteration {i_iter} ; u = {u_boat[i_iter]}')
    write_u_file(filename=u_file, u_x=u_boat[i_iter], verbose=False)
    if u_boat[i_iter] > 0.001:
        ftt = boat_length/u_boat[i_iter]
    else:
        print('Warning, very slow boat')
        ftt = 0.1
    a = int(np.ceil(M0_ftt_length*ftt))
    b = int(np.ceil((M1_ftt_length+i_iter/2)*ftt))
    c = int(np.ceil(M1_discarded_ftt_length*ftt))
    u_app = compute_app_speed(np.array([-u_boat[i_iter], 0]), u_wind)
    # print(f"endtime = {a} and {b}")
    edit_controlDict(controlDict_0_file, endTime=a)
    edit_controlDict(controlDict_1_file, endTime=b)
    lauch_computation(opf_script)
    D = compute_average_drag(drag_file, skip=c, plot_res=False)
    D_all[i_iter] = D
    wind_force = compute_wind_force(u_app, verbose=False, plot_res=False)
    T = -wind_force[0]
    T_all[i_iter] = T
    print(f'---> u_app = {u_app} ; ')
    print(f'---> Wind thrust T = {T} ; water drad D = {D}')
    print(f'---> Iteration completed in {time.time()-t1} s')
    relax_factor = relax_factors[i_iter]
    u_boat[i_iter+1] = relax_factor*u_boat[i_iter]*-T/D+(1-relax_factor)*u_boat[i_iter]
    t1 = time.time()
    if np.isnan(u_boat[i_iter+1]):
        print(f'NaN encountered. Exiting')
        print(f"u : {u_boat}")
    if (abs(u_boat[i_iter+1]-u_boat[i_iter]))/u_boat[i_iter] < threshold:
        print('---> Convergence reached. Exiting the loop')
        converged = True
        break

print(f"{i_iter+1} iterations done in {time.time()-t0} ")


# Result export and plotting 

filename = 'iterative_process_result'
if converged:
    description = f'iter u_b, T, D\n# Process converged in {i_iter-1} iterations\n# u_w={u_wind}'
else:
    description = f'iter u_b, T, D\n# Process did not converge\n"# u_w={u_wind}'
matrix = np.concatenate([np.reshape(np.arange(max_iter+1), [max_iter+1, 1]), np.reshape(u_boat, [max_iter+1, 1]), np.reshape(T_all, [max_iter+1, 1]), np.reshape(D_all, [max_iter+1, 1])], axis=1)
export_latex_curves(filename, matrix, description, verbose=True)

filename = '650_m1/postProcessing/rigidBodyMotionDisplacement/tq'
undersample_file(filename)
filename = '650_m1/postProcessing/forces/0/force.dat'
undersample_file(filename)


fig, ax=plt.subplots(1, dpi=150)
ax.plot(np.arange(len(u_boat)), u_boat)
ax.set_title('velocity boundary condition')
ax.set_ylabel('Velocity amplitude')
ax.set_xlabel('Iteration')

plt.show()

# rho = 1000
# Sp = 0.9
# Spp = 8.64

# u = 5.12
# Dp = 1460 # col 5 
# Dpp = 375 # col 8

# Cdp = 2*Dp/(rho * u**2 * Sp)
# Cdpp = 2*Dpp/(rho * u**2 * Spp)

# print(f'u = {u}, Cdp = {Cdp}, Cdpp = {Cdpp}')


# thetas = np.linspace(30,150, 10)*np.pi/180
# us = np.zeros_like(thetas)
# k = 0
# for theta in thetas:
#     u_wind = np.array([np.cos(theta), np.sin(theta)])*10
#     us[k] = estimate_initial_velocity(u_wind, verbose=False, plot_res=False)
#     # print(f'theta = {theta}, u = {u}')
#     k+=1

# for k in np.arange(len(thetas)):
#     print(f'({thetas[k]*180/np.pi},{us[k]})')