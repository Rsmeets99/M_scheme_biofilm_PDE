# Author: Robin Smeets 
# Email: robinsmeets99@gmail.com / r.k.h.smeets@uva.nl
# Institute: Korteweg-de Vries Institute for Mathematics - University of Amsterdam

'''
Python script containing all the tools for running the examples, e.g. creating the mesh and the M-scheme.
Code to solve the PDE: partial_t u = Delta (u^m) + beta * u
'''

import argparse
from dataclasses import dataclass
from inspect import getsourcefile
import os
from os import path as op
import pickle
import shutil

from dolfinx import io                  # For exporting data as xdmf files
from dolfinx import mesh                # For creating mesh
from dolfinx import fem                 # For setting up the FEM
from dolfinx.fem import FunctionSpace   # For setting up function spaces
import dolfinx.fem.petsc                # For setting up solvers
from mpi4py import MPI                  # For multicore processing
import numpy as np                      # Used in some functions
from petsc4py import PETSc              # For solving the systems of equations
import ufl                              # For defining functions and variational form
import ufl.operators as operators       # For using the max_value operator in the definition of L^i_n

file_name = op.abspath(getsourcefile(lambda:0))
file_dir = op.dirname(file_name)

## General data classes.
@dataclass
class geometry_param:
    dim: int        # Dimension of the domain
    h: float        # Mesh size
    x_min: float    # Left boundary interval/square
    x_max: float    # Right boundary interval/square

@dataclass
class solution_param:
    example_name: str       # Name of the example

    start_time: float       # Starting time of solution
    final_time: float       # Stopping time of solution
    dt: float               # Time-step size
    gamma: float            # Parameter as introduced in M-scheme
    M_par: float            # Parameter as intoruced in M-scheme

    m: float                # m parameter inside PME equation
    beta: float             # Beta parameter inside PME equation
    C: float                # Parameter for initial exact Barenblatt solution

    stop_crit: float        # Error tolerance between two iterations to decide on covergence
    allowed_iter: int       # Maximum amount of allowed iterations before deciding no convergence
    full_iter: bool         # If True, forces to make the maximum allowed of iterations

@dataclass
class general_param:
    save_plot: bool             # If true, saves .bp files (VTK format) for plotting of the solution
    save_results: bool          # If true, saves the data (e.g. amount of iterations and error each iteration)
    save_time_error: bool       # If true, calculates the error with the exact solution as given in equation 3.5 (v_n = 0)
    windows_dir: str  = None    # Name of directory in Windows where the Simulation folder can be situated

@dataclass
class experiment_param:
    dt_list: list[float]        # List of different dt one wants to run the example for
    h_list: list[float] = None  # List of different h one wants to run the example for
    M_list: list[float] = None  # List of different M_par one wants to run the example for


# Functions for the weak form and L^i_n
def Phi(u, m: float):
    return u**m 
    
def Phi_prime(u, m: float):
    return m*u**(m-1) 

# Functions for defining the exact Barenblatt solution
def alpha_BB(dim: float, m: float):
    return dim/(dim*(m-1) + 2)
    
def beta_BB(dim: float, m: float):
    return 1/(dim*(m-1) + 2)
    
def kappa_BB(dim: float, m: float):
    return (m-1)/(2*m) * 1/(dim*(m-1) + 2)

# Geometry class and Soltuion class
class geometry_class():
    """
    Class for generating the domain, function spaces and boundary conditions.
    """
    def __init__(self, param: geometry_param) -> None:
        
        self.param = param
        
        self.create_domain()
        self.create_spaces()

    def create_domain(self) -> None:
        n = int((self.param.x_max - self.param.x_min)/self.param.h)

        if self.param.dim == 1:
            self.domain = mesh.create_interval(MPI.COMM_WORLD, n, points = (self.param.x_min, self.param.x_max))
        if self.param.dim == 2:
            self.domain = mesh.create_rectangle(MPI.COMM_WORLD, [np.array([self.param.x_min, self.param.x_min]), np.array([self.param.x_max, self.param.x_max])], [n, n], mesh.CellType.triangle)

        self.x = ufl.SpatialCoordinate(self.domain)

    def create_spaces(self) -> None:
        el1 = ufl.FiniteElement("CG", self.domain.ufl_cell(), 1)
        el2 = ufl.FiniteElement("DG", self.domain.ufl_cell(), 1)
        mel = ufl.MixedElement([el1, el2])

        self.Z = FunctionSpace(self.domain, mel)
        self.V = FunctionSpace(self.domain, ("CG", 1)) 
        self.U = FunctionSpace(self.domain, ("DG", 1))


class solution_class():
    """
    Class containing the solution and solving for u_n and w_n through the M-scheme.
    """
    def __init__(self, geometry: geometry_class, param: solution_param) -> None:
        
        self.geometry: geometry_class = geometry
        self.param: solution_param = param

        self.alpha_BB: float = alpha_BB(self.geometry.param.dim, self.param.m)
        self.beta_BB: float = beta_BB(self.geometry.param.dim, self.param.m)
        self.kappa_BB: float = kappa_BB(self.geometry.param.dim, self.param.m)

        self.t: float = param.start_time
        
        self.i: int = 0
        self.j: int = 0
        self.k: int = 0

        self.Counter_list : list[float] = []
        self.Error_L2 : list[float] = []
        self.Error_grad : list[float] = []
        self.Error_main : list[float] = []
        self.Error_time: float = 0

        self.total_steps : int = round( (self.param.final_time - self.param.start_time) / self.param.dt)

        self.print_non_convergence : bool = False

        self.create_functions()
        self.initialize_sol()

    def create_functions(self) -> None:
        # Test and trial spaces for solving.
        self.v_trial, self.v_test = ufl.TrialFunction(self.geometry.V), ufl.TestFunction(self.geometry.V)
        self.w_trial, self.u_trial = ufl.TrialFunctions(self.geometry.Z)
        self.w_test, self.u_test = ufl.TestFunctions(self.geometry.Z)
        self.Z_sol = fem.Function(self.geometry.Z)

        # Defining all the functions for solving the system.
        self.u_n = fem.Function(self.geometry.U)
        self.u_n.name = "u_n"
        self.u_n_i = fem.Function(self.geometry.U)
        self.u_n_i.name = "u_n_i"
        self.w_n = fem.Function(self.geometry.V)
        self.w_n.name = "w_n"
        self.w_n_i = fem.Function(self.geometry.V)
        self.w_n_i.name = "w_n_i"

        # Defining functions for the exact solutions.
        self.u_exact = fem.Function(self.geometry.U)
        self.u_exact.name = "u_exact"
        self.w_exact = fem.Function(self.geometry.V)
        self.w_exact.name = "w_exact"

    def barenblatt(self, x, t):
        if self.geometry.param.dim == 1:
            return np.full(x.shape[1], 1/(t**self.alpha_BB) 
                           * np.maximum(0, self.param.C - self.kappa_BB * (np.sqrt(x[0]**2)/(t**self.beta_BB))**2 )**(1/(self.param.m - 1)))
        if self.geometry.param.dim == 2:
            return np.full(x.shape[1], 1/(self.t**self.alpha_BB) 
                           * np.maximum(0, self.param.C - self.kappa_BB * (np.sqrt(x[0]**2 + x[1]**2)/(t**self.beta_BB))**2 )**(1/(self.param.m - 1)))
    
    # The exact solution for the modified PME (expressed through the Barenblatt solution)
    def exact_sol_u(self, x):
        t = (1/(self.param.beta * (self.param.m - 1))) * np.e**(self.param.beta * (self.param.m - 1) * self.t)
        return np.e**(self.param.beta * self.t) * self.barenblatt(x,t)
    
    def exact_sol_w(self, x):
        return Phi(self.exact_sol_u(x), self.param.m)

    def initialize_sol(self) -> None: # Initializes u_n and u_n_i
        self.u_n.interpolate(self.exact_sol_u)
        self.u_n_i.interpolate(self.exact_sol_u)
        self.w_n.interpolate(self.exact_sol_w)      # Not necessary for solving system but for calculating Error_time
        self.w_n_i.interpolate(self.exact_sol_w)    # Not necessary for solving system but for calculating Error_time

    def L_func(self, u): # Defines L^i_n
        return operators.max_value(Phi_prime(u,self.param.m) + self.param.M_par * self.param.dt**self.param.gamma,
                                   2*self.param.M_par*self.param.dt**self.param.gamma)

    def bilinear_functional_uw(self): # Defines the bilinear functional for solving for u_n^i and w_n^i
        return ((1- self.param.dt * self.param.beta) * self.u_trial * self.w_test 
                + self.param.dt * ufl.inner(ufl.grad(self.w_trial), ufl.grad(self.w_test)) 
                + self.L_func(self.u_n_i) * self.u_trial * self.u_test - self.w_trial * self.u_test)
    
    def linear_functional_uw(self): # Defines the linear functional for solving for u_n^i and w_n^i
        return (self.u_n * self.w_test + self.L_func(self.u_n_i) * self.u_n_i * self.u_test - Phi(self.u_n_i, self.param.m) * self.u_test)
    
    def update_uw(self) -> None: # Executes the M-scheme to get to the next iteration of u_n^i and w_n^i
        
        Error = 1
        counter = 0
        
        Error_L2_temp = []
        Error_grad_temp = []
        Error_main_temp = []
        
        if self.param.full_iter:
            self.param.stop_crit = -1
        
        while Error > self.param.stop_crit:
            counter += 1 
            if counter > self.param.allowed_iter:
                if self.print_non_convergence:
                    print('did not converge')
                    print(f'H1 error is: {error_H1:.2E}') 
                    print(f'L2 error is: {error_L2:.2E}')
                    print(f'Semi-H1 error is: {error_L2grad:.2E}')
                counter = -1 # Signifies that the method did not converge at all.
                break

            # Solve the system        
            a_bilinear_uw = self.bilinear_functional_uw() * ufl.dx
            L_linear_uw = self.linear_functional_uw() * ufl.dx
            
            bilinear_form_uw = fem.form(a_bilinear_uw)
            linear_form_uw = fem.form(L_linear_uw)
            
            b_uw = dolfinx.fem.petsc.assemble_vector(linear_form_uw)
            b_uw.assemble()
            A_uw = dolfinx.fem.petsc.assemble_matrix(bilinear_form_uw)
            A_uw.assemble()
            
            solver = PETSc.KSP().create(self.geometry.domain.comm)
            solver.setType(PETSc.KSP.Type.PREONLY) # PREONLY means you only use the preconditioner
            solver.getPC().setType(PETSc.PC.Type.LU) # The LU preconditioner just solves the system using LU factorization.
            solver.setOperators(A_uw)
            
            solver.solve(b_uw, self.Z_sol.vector)
            
            w_split, u_split = self.Z_sol.sub(0).collapse(), self.Z_sol.sub(1).collapse()
            
            u_split.x.array[:] = np.maximum(np.zeros(len(u_split.x.array)), u_split.x.array)

            # Calculate the new error
            L2_u_error = fem.form(ufl.inner(self.L_func(self.u_n_i)*(u_split - self.u_n_i), u_split - self.u_n_i) * ufl.dx)
            L2grad_w_error = fem.form(ufl.inner(ufl.grad(w_split - self.w_n_i), ufl.grad(w_split - self.w_n_i)) * ufl.dx)

            error_localL2 = fem.assemble_scalar(L2_u_error)
            error_localL2grad = fem.assemble_scalar(L2grad_w_error)
            
            error_L2 = np.sqrt(self.geometry.domain.comm.allreduce(error_localL2, op=MPI.SUM))
            error_L2grad = np.sqrt(self.geometry.domain.comm.allreduce(error_localL2grad, op=MPI.SUM))
            
            Error_L2_temp.append(error_L2)
            Error_grad_temp.append(error_L2grad)
            
            error_H1 = np.sqrt(error_L2**2 + self.param.dt*error_L2grad**2)
            Error_main_temp.append(error_H1)
            
            Error = error_H1
            
            # Set new u_n^i and w_n^i
            self.u_n_i.x.array[:] = u_split.x.array
            self.w_n_i.x.array[:] = w_split.x.array

        # Set new u_n and w_n
        self.u_n.x.array[:] = self.u_n_i.x.array
        self.w_n.x.array[:] = self.w_n_i.x.array

        self.Error_L2.append(Error_L2_temp)
        self.Error_grad.append(Error_grad_temp)
        self.Error_main.append(Error_main_temp)
        self.Counter_list.append(counter)

    def time_integration(self) -> float:
        '''Compute the time integrals in equation 3.5 using the Simpson 3/8 quadrature rule'''
        int_error: float = 0

        step_sizes: list[float] = [i * self.param.dt / 3 for i in range(0,4)]
        quadrature_weights: list[float] = [1, 3, 3, 1]

        for weight, i in enumerate(quadrature_weights):

            self.t = self.t + step_sizes[i]

            self.u_exact.interpolate(self.exact_sol_u)
            self.w_exact.interpolate(self.exact_sol_w)

            L2_u_error_true = fem.form(ufl.inner(self.u_exact - self.u_n, self.u_exact - self.u_n) * ufl.dx)
            L2_w_error_true = fem.form(ufl.inner(self.w_exact - self.w_n, self.w_exact - self.w_n) * ufl.dx)

            error_localL2_u_true = fem.assemble_scalar(L2_u_error_true)
            error_localL2_w_true = fem.assemble_scalar(L2_w_error_true)
            
            error_L2_u_true = self.geometry.domain.comm.allreduce(error_localL2_u_true, op=MPI.SUM)
            error_L2_w_true = self.geometry.domain.comm.allreduce(error_localL2_w_true, op=MPI.SUM)

            int_error += weight *  (error_L2_u_true + error_L2_w_true)

        return (self.param.dt / 8) * int_error

    def calculate_error_time(self, n: int) -> None:
        self.t = (n - 1)  * self.param.dt + self.param.start_time # Exact solution is integrated from t_{n-1} to t_{n}
        self.Error_time = self.Error_time + self.time_integration()

    def plot(self, n: int) -> None: # Save the solutions to .bp file in VTX format for later plotting in e.g. ParaView

        PATH = op.join(file_dir, 'Simulations', self.param.example_name)
        if not op.exists(PATH):
            os.makedirs(PATH)
        
        self.t = n * self.param.dt + self.param.start_time
        file_name_plot_u = op.join(PATH, f'simulation_u.bp')

        if n == 0: 
            self.vtx_writer_u = io.VTXWriter(self.geometry.domain.comm, file_name_plot_u, [self.u_n], engine="BP4")
        
        self.vtx_writer_u.write(self.t)

        if n == self.total_steps:
            self.vtx_writer_u.close()

    def plot_exact(self, n: int) -> None:
        PATH = op.join(file_dir, 'Simulations', self.param.example_name)
        if not op.exists(PATH):
            os.makedirs(PATH)
        
        self.t = n * self.param.dt + self.param.start_time
        file_name_plot_u = op.join(PATH, f'simulation_u_exact.bp')

        self.u_exact.interpolate(self.exact_sol_u)

        if n == 0: 
            self.vtx_writer_u_exact = io.VTXWriter(self.geometry.domain.comm, file_name_plot_u, [self.u_exact], engine="BP4")
        
        self.vtx_writer_u_exact.write(self.t)

        if n == self.total_steps:
            self.vtx_writer_u_exact.close()

    def save_results(self) -> None: # Save the results in a pickle file
        PATH = op.join(file_dir, 'Results', self.param.example_name)
        if not op.exists(PATH):
            os.makedirs(PATH)
        file_name_results = op.join(PATH, f'results_{self.i}_{self.j}_{self.k}.pickle')
        with open(file_name_results, 'wb') as file:
            pickle.dump(self.Counter_list, file)
            pickle.dump(self.Error_L2, file)
            pickle.dump(self.Error_grad, file)
            pickle.dump(self.Error_main, file)
            pickle.dump(self.Error_time, file)

    def move_to_windows(self, windows_dir: str) -> None: # Copy the results to a given windows_directory from within WSL
        
        local_dir = op.join(file_dir, 'Simulations', self.param.example_name)
        target_dir = op.join(windows_dir, self.param.example_name)

        if not op.exists(target_dir):
            print('Created new target directory')
            os.makedirs(target_dir)

        file_name_u = op.join(local_dir, f'simulation_u.bp')
        target_file_name_u = op.join(target_dir, f'simulation_u.bp')

        shutil.copytree(file_name_u, target_file_name_u, dirs_exist_ok=True)

        file_name_u = op.join(local_dir, f'simulation_u_exact.bp')
        target_file_name_u = op.join(target_dir, f'simulation_u_exact.bp')

        shutil.copytree(file_name_u, target_file_name_u, dirs_exist_ok=True)


## Complete numerical scheme
def pme_M_scheme(solution: solution_class, 
                 parameters: general_param, 
                 parsed_args: argparse.Namespace,
                 print_progress: bool = False) -> None:
    """
    Execute the entire numerical scheme for solving the pme model

    Parameters
    ----------
    solution : solution_class
        Solution class that contains the solution and functions for solving for new time steps
    parameters : general_param
        General class containing information on whether or not to save data and/or save the solutions for plotting
    """

    if parameters.save_plot:
        solution.plot(0)
        solution.plot_exact(0)

    # Loop over time steps.
    for n in range(solution.total_steps):

        # Solve for (u_n, w_n) using M-scheme.
        solution.update_uw()

        if solution.Counter_list[-1] == -1:
            if print_progress:
                print('Stopped example as did not converge')
            if parameters.save_plot:
                solution.plot(solution.total_steps) # Make sure that writer is closed.
                solution.plot_exact(solution.total_steps) # Make sure that writer is closed.
            break

        if parameters.save_plot:
            solution.plot(n+1)
            solution.plot_exact(n+1)

        if parameters.save_time_error:
            solution.calculate_error_time(n+1)

        if print_progress:
            print(f'\rScheme is at {(n+1)/solution.total_steps * 100:.2f}%.', end = '')
            
    if parameters.save_results:
        if print_progress:
            print('\nSaving results.')
        solution.save_results()

    if parsed_args.movewsl and parameters.save_plot:
        if print_progress:
            print('Copy simulation data to Windows from WSL2')
        solution.move_to_windows(parameters.windows_dir)


if __name__ == '__main__':

    import pme_config as config
    print(f"DOLFINx version: {dolfinx.__version__}")
    geometry_config = config.geometry_param_1D_avg_iter_pme
    solution_config = config.solution_param_1D_avg_iter_pme
    general_config = config.general_param_1D_avg_iter_pme 

    geometry = geometry_class(geometry_config)
    solution = solution_class(geometry, solution_config)
