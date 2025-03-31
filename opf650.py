import numpy as np 
import matplotlib.pyplot as plt 
# import sys 
import os 
import subprocess 
import fileinput 



def compute_average_drag(filename, skip=0, plot_res=False):
    forces = np.loadtxt(filename)
    t = forces[:,0]
    n_points = len(t)
    firstPoints = np.sum(t<=skip)
    drag = forces[firstPoints:,1]
    t = forces[firstPoints:,0]
    if plot_res:
        fig, ax = plt.subplots(1, dpi=100)
        ax.plot(t, drag)
        # ax.set_ylim([-2000,0])
        plt.show()
        # print(len(drag))
    return np.mean(drag)

def write_u_file(filename='U', u_x=0, u_y=0, u_z=0, verbose=True):
    if verbose:
        print(f'writing U of amplitude {np.linalg.norm([u_x, u_y, u_z])} into {filename}')
    if os.path.exists(filename):
        os.remove(filename)
    foamUFile = f"\nFoamFile\n{{\n    version     2.0;\n    format      ascii;\n    class       volVectorField;\n    object      U;\n}}\n\ndimensions      [0 1 -1 0 0 0 0];\n\ninternalField   uniform ({-u_x} {-u_y} {-u_z});\n\nboundaryField\n{{\n  #includeEtc \"caseDicts/setConstraintTypes\"\n    inlet\n    {{\n        type            fixedValue;\n        value           $internalField;\n    }}\n\n    outlet\n    {{\n        type            outletPhaseMeanVelocity;\n        alpha           alpha.water;\n        Umean           {np.linalg.norm([u_x, u_y, u_z])};\n        value           $internalField;\n    }}\n\n    atmosphere\n    {{\n        type            pressureInletOutletVelocity;\n        tangentialVelocity $internalField;\n        value           uniform (0 0 0);\n    }}\n\n    \"coque|quille|foils\"\n    {{\n        type            noSlip;\n    }}\n}}\n"

    with open(filename, 'w') as f:
        f.write(foamUFile)
        f.close()
    return True

def edit_controlDict(filename, endTime=1):
    with fileinput.FileInput(filename, inplace = True, backup ='.bak') as f: 
        for line in f: 
            if "endTime" in line[0:10]: 
                print(f"endTime         {endTime};", end ='\n') 
            else: 
                print(line, end ='')

def lauch_computation(opf_script):
    subprocess.run(opf_script)
    return True


def compute_app_speed(u_boat, u_wind):
    u_app = u_wind-u_boat
    return u_app



    # sail angle, camber en deg par rapport  l'axe du bateau. cf figure du rapport. 
    # en upwind, la voile marche comme un profil aéro
    # en downwind, c'est le drag qui pousse. Pour savoir dans quelle config on est, on calcule les deux et on regarde la plus grande. 


def compute_wind_force(u_app, sail_angle, verbose=True, plot_res=False):
    rho=1.2
    S = 34 # en m2
    camber = 20
    # sail_angle = sail_angle + camber
    sail_angle *= np.pi/180
    camber *= np.pi/180
    sail = np.array([np.cos(sail_angle), np.sin(sail_angle)])
    rotation_mat = np.array([[np.cos(sail_angle), -np.sin(sail_angle)], [np.sin(sail_angle), np.cos(sail_angle)]])
    print(f"rotation mat = {rotation_mat}")
    # upwind hypothesis 
    Cl = 0.8   # réglage au près - upwind 
    Cd = 0.02  # réglage au près - upwind
    u = u_app@sail.T
    Fl_upwind = 0.5*rho*u**2*S*Cl
    Fd_upwind = 0.5*rho*u**2*S*Cd
    F_upwind = rotation_mat@np.array([Fd_upwind, Fl_upwind]).T
    F_upwind_pas_rot = np.array([Fd_upwind, Fl_upwind])
    if verbose:
        print(f'---> Upwind mode, u_app = {u_app}, u*sail = {u}, F = {F_upwind}')
    # pas upwind hypothesis 
    bullshit_coeff = 0.5
    sail = np.array([-np.sin(sail_angle), np.cos(sail_angle)])
    u = u_app@sail.T
    stagnation_pressure = 0.5*rho*u**2*bullshit_coeff 
    Fl_downwind = 0
    Fd_downwind = S * stagnation_pressure  
    F_downwind = rotation_mat@np.array([Fl_downwind, Fd_downwind]).T
    if verbose:
        print(f'---> Downwind mode, u_app = {u_app}, u*sail = {u}, F = {F_downwind}')
    F = F_upwind
    if abs(F_downwind[0]) > abs(F_upwind[0]):
        F = F_downwind
    if verbose:
        print(f"F = {F}")
    if F[0] > 0:
        print(f"WARNING The wind thrust is not pushing forward")
    if plot_res:
        fig, ax = plt.subplots(1, dpi=100)
        ax.set_xlim(-1, 1)
        ax.set_ylim(-0.5, 0.5)
        ax.set_aspect('equal')
        # Bateau, voile
        boat_length = 1.0
        boat_width = 0.3
        boat_shape = np.array([[-boat_length / 2, 0], [0, -boat_width / 2], [boat_length / 2, -boat_width / 2], [boat_length / 2, boat_width / 2], [0, boat_width / 2], [-boat_length / 2, 0]])
        ax.add_patch(plt.Polygon(boat_shape, color='gray', alpha=0.6))
        sail_x = [0, np.cos(sail_angle) * 0.5]
        sail_y = [0, np.sin(sail_angle) * 0.5]
        ax.plot(sail_x, sail_y, color='black', linewidth=4)
        # Normalisation pour la mise à l'échelle des flèches
        norm_fact = boat_length / max(np.linalg.norm(F_upwind), np.linalg.norm(F_downwind))*0.8
        # Vent apparent
        u_app = u_app/np.linalg.norm(u_app)*0.2
        ax.arrow(-0.75, 0, u_app[0], u_app[1], head_width=0.05, color='black')
        ax.arrow(-0.75, 0.1, u_app[0], u_app[1], head_width=0.05, color='black')
        ax.arrow(-0.75, -0.1, u_app[0], u_app[1], head_width=0.05, color='black')
        # Forces
        ax.arrow(0, 0, F_upwind[0] * norm_fact, F_upwind[1] * norm_fact, head_width=0.03, color='red', linestyle='solid', label='upwind mode')
        ax.arrow(0, 0, F_upwind_pas_rot[0] * norm_fact, F_upwind_pas_rot[1] * norm_fact, head_width=0.03, color='red', linestyle='dotted', label='upwind mode')
        ax.arrow(0, 0, F_downwind[0] * norm_fact, F_downwind[1] * norm_fact, head_width=0.03, color='blue', linestyle='solid', label='downwind mode')
        ax.legend()
        ax.set_xlabel("x ")
        ax.set_ylabel("y ")
        ax.set_title("Forces et vent sur le bateau")
        plt.grid()
        plt.show()
    return F

# compute_wind_force(u_app, sail_angle, verbose=verbose, plot_res=plot_res)

def remove_alternate_lines(file_path):
    # Lire toutes les lignes du fichier
    with open(file_path, 'r', encoding='utf-8') as file:
        lines = file.readlines()
    
    # Garder une ligne sur deux
    filtered_lines = lines[::2]
    
    # Réécrire le fichier avec les lignes filtrées
    with open(file_path, 'w', encoding='utf-8') as file:
        file.writelines(filtered_lines)
        
def undersample_file(file_path, N=200):
    # Lire toutes les lignes du fichier
    with open(file_path, 'r', encoding='utf-8') as file:
        lines = file.readlines()
    
    # Vérifier si le fichier est plus petit que N lignes
    if len(lines) <= N:
        return
    
    # Sélectionner N lignes réparties uniformément
    indices = np.linspace(0, len(lines) - 1, N, dtype=int)
    filtered_lines = [lines[i] for i in indices]
    
    # Réécrire le fichier avec les lignes filtrées
    with open(file_path, 'w', encoding='utf-8') as file:
        file.writelines(filtered_lines)
