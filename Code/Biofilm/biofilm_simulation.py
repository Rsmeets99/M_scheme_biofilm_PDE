import time

## Import support script
time_package0 = time.time()
from biofilm_utils import *
time_package1 = time.time()
print(f'Importing support script took {time_package1-time_package0:.1f} seconds')

time_setup0 = time.time()

# Options for plotting (will make the corresponding .xdmf and .h5 files for use in ParaView)
plotting = True # Plots numerical solution
print_progress = True # Prints max and min values at each time step, as well as amount of iterations it took

# Convergence criteria M-scheme
iter_M = 100 # Maximum amount of iterations before saying it did not converge
stop_crit = 1e-5 # Tolerance for error to say it has converged

# General parameters
t_start = 0
T = 1.5
dt = 10**(-2.5)
h = 1/200 

gamma = 1/3
M_par = 1e-3

# Parameters for initial conditions
# Note: make sure that you choose your initial blobs in such a way that they won't overlap with as a result that the initial condition goes above 1.
x0, y0 = 0.0, 0.0
x1, y1 = -0.3, -0.3
x2, y2 = 0.3, 0.3
x3, y3 = -0.3, -0.3
x4, y4 = -0.6, -0.3
xlist = [x1,x2] # x coordinates centers initial blobs
ylist = [y1,y2] # y coordinates centers initial blobs, if dim = 1 you can ignore this one, as long as it still exists (even if empty)
r = 0.2
height = 0.95

# Deciding between PDE-PDE and PDE-ODE case
mu = 1

# Defining the domain
dim = 1
if dim == 1:
    sizes = -1, 1 # x_min and x_max, for 2D use x_min, x_max, y_min, y_max
if dim == 2:
    sizes = -1, 1, -1, 1

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

# Generating mesh for domain
domain, x = create_domain(h, dim, sizes, use_gmsh= False)

# Creating solution spaces, and trial/test functions
Z, V, U, w_trial, u_trial, w_test, u_test, Z_sol, v_trial, v_test = create_space_sol(domain)

# Creating boundary condition for v_n in the PDE-PDE case
bc2 = boundary_conditions_v(dim, U, V, domain, x, sizes)

# Creating the classes for efficiency
L_func, expr_u = create_classes(dim, r, height, M_par, gamma, dt, un_max_estimate_4, L_estimate_4, delta_1)

# Creating functions
Diffusion_v, u_n, u_n_i, w_n, w_n_i, v_n, uh, wh, vh, v_proj = create_functions(domain, U, V, delta_2)

time_setup1 = time.time()
print(f'Set up took {time_setup1-time_setup0:.1f} seconds')

time_running0 = time.time()
Run_M_scheme(M_par = M_par, dt = dt, h = h, T = T, t_start = t_start, gamma = gamma, dim = dim, sizes = sizes, r = r, height = height, xlist = xlist, ylist = ylist, mu = mu, un_max_estimate_4 = un_max_estimate_4, L_estimate_4 = L_estimate_4, Phi_un_max_4 = Phi_un_max_4, delta_1 = delta_1, kappa_2 = kappa_2, kappa_3 = kappa_3, kappa_4 = kappa_4, kappa_1 = kappa_1, domain = domain, U = U, u_trial = u_trial, u_test = u_test, w_trial = w_trial, w_test = w_test, v_trial = v_trial, v_test = v_test, v_proj = v_proj, u_n_i = u_n_i, w_n_i = w_n_i, v_n = v_n, u_n = u_n, w_n = w_n, uh = uh, wh = wh, vh = vh, Z_sol = Z_sol, Diffusion_v = Diffusion_v, bc2 = bc2, L_func= L_func, expr_u= expr_u, iter_M = iter_M, stop_crit = stop_crit, plotting = True, full_iter = False, save_data = False, data_return = False, pickle_file = False, print_progress = print_progress)
time_running1 = time.time()
print(f'running took {time_running1-time_running0:.1f} seconds')
print(f'Total took {time_running1-time_package0:.1f} seconds')