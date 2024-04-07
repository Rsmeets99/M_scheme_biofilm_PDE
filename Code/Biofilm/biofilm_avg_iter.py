import time
import matplotlib.pyplot as plt

## Import support script
from biofilm_utils import *

# Convergence criteria M-scheme
iter_M = 400 # Maximum amount of iterations before saying it did not converge
stop_crit = 1e-5 # Tolerance for error to say it has converged

# General parameters
t_start = 0
T = 1.2
dt = 10**(-1)

gamma = 1/4
 
# Parameters h and M for h dependency
h_list = [1/i for i in range(10,210,10)]
hinv_list = [i for i in range(10,210,10)]

M_list = [1e-7, 5e-5, 5e-4, 5e-3]

# Deciding between PDE-PDE and PDE-ODE case
mu = 0

# Defining size domain
dim = 1
if dim == 1:
    sizes = -2, 2
if dim == 2:
    sizes = -2, 2, -2, 2

# Parameters for initial conditions
x0, y0 = -0.4, 0.5
x1, y1 = -0.3, -0.3
x2, y2 = 0.3, 0.3
x3, y3 = -0.3, -0.3
x4, y4 = -0.6, -0.3
xlist = [x1,x2]
ylist = [y1,y2]
r = 0.2
height = 0.9

# Defining size domain
dim = 1
sizes = -2, 2 # x_min and x_max, for 2D use x_min, x_max, y_min, y_max
# sizes = -2, 2, -2, 2

# Defining function constants
kappa_1 = 0.4
kappa_2 = 0.01
kappa_3 = 1
kappa_4 = 0.42
delta_1 = 1*10**(-6)
delta_2 = 0.2

# Take un_max_estimate_4 big enough. Can estimate it using proven proposition
# TODO add script to calculate un_max_estimate?
un_max_estimate_4 = 0.993
L_estimate_4 = Phi_prime_4(un_max_estimate_4, delta_1)
Phi_un_max_4 = Phi_4_np(un_max_estimate_4, delta_1)

def generate_data(h_list, M_list):
    data_avg = np.zeros((len(M_list),len(h_list)))

    time_hstart = time.time()
    for i in range(len(h_list)):
        time_h0 = time.time()
        # Generating mesh for domain
        h = h_list[i] # mesh size
        domain, x = create_domain(h, dim, sizes)

        # Creating solution spaces, and trial/test functions
        Z, V, U, w_trial, u_trial, w_test, u_test, Z_sol, v_trial, v_test = create_space_sol(domain)

        # Creating boundary condition for v_n in the PDE-PDE case
        bc2 = boundary_conditions_v(dim, U, V, domain, x, sizes)
        
        # Creating functions
        Diffusion_v, u_n, u_n_i, w_n, w_n_i, v_n, uh, wh, vh, v_proj = create_functions(domain, U, V, delta_2)
        
        for j in range(len(M_list)):
            # Creating the classes for efficiency
            L_func, expr_u = create_classes(dim, r, height, M_list[j], gamma, dt, un_max_estimate_4, L_estimate_4, delta_1)
            
            EL2_temp, EL2grad_temp, EH1_temp, Counter_array, avg_iter = Run_M_scheme(M_par = M_list[j], dt = dt, h = h, T = T, t_start = t_start, gamma = gamma, dim = dim, sizes = sizes, r = r, height = height, xlist = xlist, ylist = ylist, mu = mu, un_max_estimate_4 = un_max_estimate_4, L_estimate_4 = L_estimate_4, Phi_un_max_4 = Phi_un_max_4, delta_1 = delta_1, kappa_2 = kappa_2, kappa_3 = kappa_3, kappa_4 = kappa_4, kappa_1 = kappa_1, domain = domain, U = U, u_trial = u_trial, u_test = u_test, w_trial = w_trial, w_test = w_test, v_trial = v_trial, v_test = v_test, v_proj = v_proj, u_n_i = u_n_i, w_n_i = w_n_i, v_n = v_n, u_n = u_n, w_n = w_n, uh = uh, wh = wh, vh = vh, Z_sol = Z_sol, Diffusion_v = Diffusion_v, bc2 = bc2, L_func= L_func, expr_u= expr_u, iter_M = iter_M, stop_crit = stop_crit, plotting = False, full_iter = False, save_data = False, data_return = True, pickle_file = False, print_progress = False)

            data_avg[j,i] = avg_iter

        time_h1 = time.time()
        print(f'Done with 1/h = {int(1/h):.0f}. It took {time_h1 - time_h0:.0f} seconds. Running for a total of {time_h1 - time_hstart:.0f} seconds now')
    return data_avg

def plot_hdep_avg(data_avg, hinv_list, M_list, dt, dim, mu, T):
    for j in range(len(M_list)):
        plt.plot(hinv_list, data_avg[j], '--.', label = f'M = {M_list[j]:.1E}')

    plt.xlabel('Mesh size 1/h')
    plt.ylabel('Average amount of iterations')
    plt.title(f'$\\tau$ = {dt:.1E} and T = {T:.1f}')
    plt.grid()
    plt.legend()
    plt.savefig(f'biofilm_hdep_avg_dim{dim:.0f}D_mu{mu:.0f}_dt{dt:.2f}_gamma{gamma:.2f}.png')
    plt.clf()
    return 0


def main_nodata(hinv_list, M_list, h_list, dt, dim, mu, save_data = True):
    data_avg = generate_data(h_list, M_list)
    if save_data == True:
        with open(f'biofilm_hdep_{dim:.0f}D_mu{mu:.0f}_dt{dt:.2f}.pickle', 'wb') as file:
            pickle.dump(data_avg, file)
            pickle.dump([M_list, hinv_list, M_list, dt, dim, mu, T], file)

    plot_hdep_avg(data_avg, hinv_list, M_list, dt, dim, mu, T)
    return 0
    
def main_withdata(dim,dt):
    data = []
    with (open(f'biofilm_hdep_{dim:.0f}D_mu{mu:.0f}_dt{dt:.2f}.pickle', 'rb')) as openfile:
        while True:
            try:
                data.append(pickle.load(openfile))
            except EOFError:
                break
    data_avg, data_gen = data
    M_list, h_list, hinv_list, dt, dim, T, t_start = data_gen

    plot_hdep_avg(data_avg, hinv_list, M_list, dt, dim, mu, T)
    return 0
    
main_nodata(hinv_list, M_list, h_list, dt, dim, mu, save_data = True)


# SMALL_SIZE = 8
# MEDIUM_SIZE = 10
# BIGGER_SIZE = 12

# plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
# plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
# plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
# plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
# plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
# plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
# plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title