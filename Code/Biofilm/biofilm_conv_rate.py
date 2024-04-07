import time
import matplotlib.pyplot as plt

## Import support script
from biofilm_utils import *

# Convergence criteria M-scheme
iter_M = 3 # We force it to make 3 iterations to compute the convergence rate
stop_crit = 1e-5 # Tolerance for error to say it has converged (we ignore this now as we force it to make 3 iterations either way)

# General parameters
t_start = 0
T = 0.5
h = 1e-4
dt_list = [10**(-i/4) for i in range(4,11)]

M_par = 5e-4
gamma = 1/4

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
x5, y5 = 0, 0
xlist = [x1,x2]
ylist = [y1,y2]
r = 0.2
height = 0.9

# Defining function constants
kappa_1 = 0.4
kappa_2 = 0.01
kappa_3 = 1
kappa_4 = 0.42
delta_1 = 1*10**(-6)
delta_2 = 0.2

# Take un_max_estimate_4 big enough. Can estimate it using proven proposition
# TODO add script to calculate un_max_estimate?
un_max_estimate_4 = 0.9999
L_estimate_4 = Phi_prime_4(un_max_estimate_4, delta_1)
Phi_un_max_4 = Phi_4_np(un_max_estimate_4, delta_1)

def geo_mean_log(iterable):
    return np.log10(np.exp(np.log(iterable).mean()))

def generate_data(dt_list):
    time_startdata = time.time()
    data = np.zeros(len(dt_list))
    # Generating mesh for domain
    domain, x = create_domain(h, dim, sizes)

    # Creating solution spaces, and trial/test functions
    Z, V, U, w_trial, u_trial, w_test, u_test, Z_sol, v_trial, v_test = create_space_sol(domain)

    # Creating boundary condition for v_n in the PDE-PDE case
    bc2 = boundary_conditions_v(dim, U, V, domain, x, sizes)
    
    # Creating functions
    Diffusion_v, u_n, u_n_i, w_n, w_n_i, v_n, uh, wh, vh, v_proj = create_functions(domain, U, V, delta_2)
    
    for i in range(len(dt_list)):
        time_h0 = time.time()

        # Creating the classes for efficiency
        L_func, expr_u = create_classes(dim, r, height, M_par, gamma, dt_list[i], un_max_estimate_4, L_estimate_4, delta_1)
    
        EL2_temp, EL2grad_temp, EH1_temp, Counter_array, avg_iter = Run_M_scheme(M_par = M_par, dt = dt_list[i], h = h, T = T, t_start = t_start, gamma = gamma, dim = dim, sizes = sizes, r = r, height = height, xlist = xlist, ylist = ylist, mu = mu, un_max_estimate_4 = un_max_estimate_4, L_estimate_4 = L_estimate_4, Phi_un_max_4 = Phi_un_max_4, delta_1 = delta_1, kappa_2 = kappa_2, kappa_3 = kappa_3, kappa_4 = kappa_4, kappa_1 = kappa_1, domain = domain, U = U, u_trial = u_trial, u_test = u_test, w_trial = w_trial, w_test = w_test, v_trial = v_trial, v_test = v_test, v_proj = v_proj, u_n_i = u_n_i, w_n_i = w_n_i, v_n = v_n, u_n = u_n, w_n = w_n, uh = uh, wh = wh, vh = vh, Z_sol = Z_sol, Diffusion_v = Diffusion_v, bc2 = bc2, L_func= L_func, expr_u= expr_u, iter_M = iter_M, stop_crit = stop_crit, plotting = False, full_iter = True, save_data = False, data_return = True, pickle_file = False, print_progress = False)

        error_data = EH1_temp[0]
        temp = [error_data[i+1]/error_data[i] for i in range(len(error_data)-1)]
        conv_rate = geo_mean_log(temp)
        data[i] = conv_rate
        time_h1 = time.time()
        print(f'We are now done with dt = {dt_list[i]:.2E}. It took {time_h1-time_h0:.0f} seconds. A total of {time_h1-time_startdata:.0f} seconds.')
    return data

dt_log_list = np.log10(np.array(dt_list))

data = generate_data(dt_list)
cutoff = 3

m,b = np.polyfit(dt_log_list[cutoff:], data[cutoff:], 1)

plt.plot(dt_log_list, data, 'r.', label = '$\log_{10}(\\alpha)$')
plt.plot(dt_log_list, m*dt_log_list + b, '--', label = f'slope = {m:.2f}')
plt.xlabel('Log of time-step dt')
plt.ylabel('Log of convergence rate $\\alpha$')
plt.title(f'h = {h:.1E}, $\\gamma$ = {gamma:.2f}, {dim:.0f}D, $\\mu$ = {mu:.0f}')
plt.grid()
plt.legend()
plt.savefig(f'biofilm_conv_rate{dim:.0f}D_mu{mu:.0f}_h{h:.1E}_gamma{gamma:.2f}_Mpar{M_par:.1E}.png')
plt.clf()

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



















