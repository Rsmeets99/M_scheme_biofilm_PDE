from dataclasses import dataclass
from inspect import getsourcefile
import os
from os import path as op
import pickle
from typing import Any

from dolfinx import io # For exporting data as xdmf files
from dolfinx import mesh # For creating mesh
from dolfinx.mesh import locate_entities, meshtags # For boundary conditions
from dolfinx import fem # For setting up the FEM
from dolfinx.fem import FunctionSpace # Function space
from mpi4py import MPI # For multicore processing
import numpy as np # Used in some functions
from petsc4py import PETSc # For solving the systems of equations
from petsc4py.PETSc import ScalarType # For defining a scalar type function
import ufl # For defining functions and variational form

file_name = op.abspath(getsourcefile(lambda:0))
file_dir = op.dirname(file_name)

## General data classes.

@dataclass
class geometry_param:
    dim: int
    h: float
    x_min: float
    x_max: float

@dataclass
class solution_param:
    example_name: str

    start_time : float
    final_time: float
    dt: float
    gamma: float
    M_par: float

    stop_crit: float
    allowed_iter: int
    full_iter: bool

    r: float
    height: float
    xlist: list[float]
    ylist: list[float]

    un_max_estimate: float


@dataclass
class general_param:
    save_plot: bool
    save_results: bool

# Functions for the weak form
def f_func(v, kappa_2, kappa_3, kappa_4):
    return kappa_3 * (v / (kappa_2 + v)) - kappa_4
    
def g_func(v, kappa_2, kappa_1):
    return -1*kappa_1 * v / (kappa_2 + v)

def Phi_4(u, delta_1): # ufl function
    return delta_1 * ( (18*u**2 - 30*u + 13)/(3*(1-u)**3) + u + 4*ufl.operators.ln(1-u) - 13/3)

def Phi_4_np(u, delta_1): # numpy function used for potentially calculating Phi_un_max_4 from un_max_estimate 
    return delta_1 * ( (18*u**2 - 30*u + 13)/(3*(1-u)**3) + u + 4*np.log(1-u) - 13/3)

def Phi_prime_4(u, delta_1): # works both for ufl and numpy
    return delta_1* u**4 / (1-u)**4

def Phi_4_reg(u, un_max_estimate_4, L_estimate_4, Phi_un_max_4, delta_1): 
    return ufl.conditional(ufl.gt(u, un_max_estimate_4), L_estimate_4 * (u-un_max_estimate_4) + Phi_un_max_4, Phi_4(u, delta_1) )
    
def Phi_prime_4_reg(u, un_max_estimate_4, L_estimate_4, delta_1): 
    return ufl.conditional(ufl.gt(u, un_max_estimate_4), L_estimate_4, Phi_prime_4(u, delta_1) )

class geometry_class():

    def __init__(self, param: geometry_param) -> None:
        
        self.param = param

        self.create_domain()
        self.create_spaces()
        self.create_boundary_condition()

    def create_domain(self) -> None:
        n = int((self.x_max - self.x_min)/self.h)
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

    def create_boundary_condition(self) -> None: 
        if self.param.dim == 1:
            dofs_V = fem.locate_dofs_geometrical(self.V, lambda x: np.isclose(x[0], self.param.x_max)) # Finds the locations close to where x = x_max (right border)
            bcV = fem.dirichletbc(ScalarType(1), dofs_V, self.V) 
        if self.param.dim == 2:
            dofs_V = fem.locate_dofs_geometrical(self.V, lambda x: np.isclose(x[1], self.param.x_max)) # Finds the locations close to where y = x_max (top border)
            bcV = fem.dirichletbc(ScalarType(1), dofs_V, self.V)
        
        self.bc = [bcV]


class solution_class():

    def __init__(self, geometry: geometry_class, param: solution_param) -> None:
        
        self.geometry: geometry_class = geometry
        self.param: solution_param = param
        
        self.i: str = 0
        self.j: str = 0
        self.k: str = 0

        self.L_max_estimate: float = Phi_prime_4(self.param.un_max_estimate, self.param.delta_1)
        self.Phi_max_estimate: float = Phi_4_np(self.param.un_max_estimate, self.param.delta_1)

        self.Counter_list : list[float] = []
        self.Error_L2 : list[float] = []
        self.Error_grad : list[float] = []
        self.Error_main : list[float] = []

        self.create_functions()
        self.initialize_sol()

        if self.param.mu == 1:
            self.set_solver_v_mu_1()

    def create_functions(self) -> None:
        # Test and trial spaces for solving.
        self.v_trial, self.v_test = ufl.TrialFunction(self.geometry.V), ufl.TestFunction(self.geometry.V)
        self.w_trial, self.u_trial = ufl.TrialFunctions(self.geometry.Z)
        self.w_test, self.u_test = ufl.TestFunctions(self.geometry.Z)
        self.Z_sol = fem.Function(self.geometry.Z)

        # Defining the function spaces for all the functions
        self.u_n = fem.Function(self.geometry.U)
        self.u_n.name = "u_n"
        self.u_n_i = fem.Function(self.geometry.U)
        self.u_n_i.name = "u_n_i"
        self.w_n = fem.Function(self.geometry.V)
        self.w_n.name = "w_n"
        self.w_n_i = fem.Function(self.geometry.V)
        self.w_n_i.name = "w_n_i"
        self.v_n = fem.Function(self.geometry.V)
        self.v_n.name = "v_n"

        # Function for the projection function used in the PDE-ODE case
        self.v_proj = fem.Function(self.geometry.V)
        self.v_proj.name = "v_proj"

        # Constant function for PDE-PDE case in weak form of v_n
        self.Diffusion_v = fem.Constant(self.geometry.domain, ScalarType(self.param.delta_2))

    def initial_u_iter(self, x, xlist, ylist): # Recursive function
        if self.geometry.param.dim == 1:
            if len(xlist) == 1:
                return np.full(x.shape[1], self.param.height*(1/self.param.r)*np.sqrt(np.maximum(0, self.param.r**2 - (x[0]-xlist[0])**2)))
            return np.full(x.shape[1], self.param.height*(1/self.param.r)*np.sqrt(np.maximum(0, self.param.r**2 - (x[0]-xlist[0])**2))) + self.initial_u_iter(x, xlist[1:], ylist)
        
        if self.geometry.param.dim == 2:
            if len(xlist)== 1:
                return np.full(x.shape[1], self.param.height*(1/self.param.r)*np.sqrt(np.maximum(0, self.param.r**2 - (x[0]-xlist[0])**2 - (x[1]-ylist[0])**2)))
            return np.full(x.shape[1], self.param.height*(1/self.param.r)*np.sqrt(np.maximum(0, self.param.r**2 - (x[0]-xlist[0])**2 - (x[1]-ylist[0])**2))) + self.initial_u_iter(x, xlist[1:], ylist[1:])

    def initialize_u(self, x) -> None:
        xlist = self.param.xlist
        ylist = self.param.ylist
        return self.initial_u_iter(x, xlist, ylist)
    
    def initialize_sol(self) -> None:
        self.u_n.interpolate(self.initialize_u)
        self.u_n_i.interpolate(self.initialize_u)
        self.v_n.x.array[:] = np.ones(len(self.v_n.x.array))

    def L_func(self, u):
        return ufl.operators.max_value(Phi_prime_4_reg(u,self.param.un_max_estimate, self.L_max_estimate, self.param.delta_1) 
                                       + self.param.M_par * self.param.dt**self.param.gamma, 
                                       2*self.param.M_par*self.param.dt**self.param.gamma)

    def bilinear_functional_uw(self):
        return ((1-f_func(self.v_n, self.param.kappa_2, self.param.kappa_3, self.param.kappa_4)*self.param.dt) * self.u_trial * self.w_test 
                + self.param.dt * ufl.inner(ufl.grad(self.w_trial), ufl.grad(self.w_test)) 
                + self.L_func(self.u_n_i) * self.u_trial * self.u_test 
                - self.w_trial * self.u_test) 
    
    def linear_functional_uw(self):
        return (self.u_n * self.w_test 
                + self.L_func(self.u_n_i) * self.u_n_i * self.u_test 
                - Phi_4_reg(self.u_n_i, self.param.un_max_estimate, self.L_max_estimate, self.Phi_max_estimate, self.param.delta_1) * self. u_test)

    def update_uw(self) -> None:

        Error = 1
        counter = 0
        
        Error_L2_temp = []
        Error_grad_temp = []
        Error_main_temp = []
        
        if self.param.full_iter:
            self.stop_crit = 0
        
        while Error > self.stop_crit:
            if counter > self.param.allowed_iter:
                print('did not converge')
                print(f'H1 error is: {error_H1:.2E}') 
                print(f'L2 error is: {error_L2:.2E}')
                print(f'Semi-H1 error is: {error_L2grad:.2E}')
                counter = -1 # Signifies that the method did not converge at all.
                break

            # Calculating new u_n^i and w_n^i
            counter += 1
                        
            a_bilinear_uw = self.bilinear_functional_uw() * ufl.dx
            L_linear_uw = self.linear_functional_uw() * ufl.dx
            
            bilinear_form_uw = fem.form(a_bilinear_uw)
            linear_form_uw = fem.form(L_linear_uw)
            
            b_uw = fem.petsc.assemble_vector(linear_form_uw)
            b_uw.assemble()
            A_uw = fem.petsc.assemble_matrix(bilinear_form_uw)
            A_uw.assemble()
            
            solver = PETSc.KSP().create(self.geometry.domain.comm)
            solver.setType(PETSc.KSP.Type.PREONLY) # PREONLY means you only use the preconditioner
            solver.getPC().setType(PETSc.PC.Type.LU) # The LU preconditioner just solves the system using LU factorization.
            solver.setOperators(A_uw)
            
            solver.solve(b_uw, self.Z_sol.vector)
            
            w_split, u_split = self.Z_sol.sub(0).collapse(), self.Z_sol.sub(1).collapse()
            
            u_split.x.array[:] = np.maximum(np.zeros(len(u_split.x.array)), u_split.x.array)

            # Can calculate error here between uh and u_n_i, or between wh and w_n_i (or both) for the tolerance
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
            
            self.w_n_i.x.array[:] = w_split.x.array
            self.u_n_i.x.array[:] = u_split.x.array 

        # Set new u_n and w_n
        self.u_n.x.array[:] = self.u_n_i.x.array
        self.w_n.x.array[:] = self.w_n_i.x.array

        self.Counter_list.append(counter)
        self.Error_L2.append(Error_L2_temp)
        self.Error_grad.append(Error_grad_temp)
        self.Error_main.append(Error_main_temp)

    def set_solver_v_mu_1(self) -> None:

        self.a_bilinear_v = (self.v_trial * self.v_test 
                             + self.Diffusion_v * self.param.dt * ufl.inner(ufl.grad(self.v_trial), ufl.grad(self.v_test))) * ufl.dx
        self.bilinear_form_v = fem.form(self.a_bilinear_v)
        
        self.L_linear_v = (self.v_n + self.param.dt*g_func(self.v_n, self.param.kappa_2, self.param.kappa_1) * self.u_n) * self.v_test * ufl.dx
        self.linear_form_v = fem.form(self.L_linear_v)
        self.b_v = fem.petsc.create_vector(self.linear_form_v)
        
        self.A_v = fem.petsc.assemble_matrix(self.bilinear_form_v, bcs= self.geometry.bc)
        self.A_v.assemble()
        
        self.solver_v = PETSc.KSP().create(self.geometry.domain.comm)
        self.solver_v.setOperators(self.A_v)
        self.solver_v.setType(PETSc.KSP.Type.PREONLY)
        self.solver_v.getPC().setType(PETSc.PC.Type.LU)

    def update_v(self) -> None:

        if self.param.mu == 0: 

            expr_v = g_func(self.v_n,self.param.kappa_2, self.param.kappa_1)*self.param.dt*self.u_n + self.v_n
            a_bilinear_proj = ufl.inner(self.v_trial, self.v_test) * ufl.dx
            L_linear_proj = ufl.inner(expr_v, self.v_test) * ufl.dx
            
            bilinear_form_proj = fem.form(a_bilinear_proj)
            linear_form_proj = fem.form(L_linear_proj)
            
            b_proj = fem.petsc.assemble_vector(linear_form_proj)
            b_proj.assemble()
            A_proj = fem.petsc.assemble_matrix(bilinear_form_proj)
            A_proj.assemble()
            
            solver = PETSc.KSP().create(self.geometry.domain.comm)
            solver.setType(PETSc.KSP.Type.PREONLY)
            solver.getPC().setType(PETSc.PC.Type.LU) 
            
            solver.setOperators(A_proj)
            
            solver.solve(b_proj, self.v_n.vector)
            self.v_n.x.scatter_forward()

            self.v_n.x.array[:] = np.maximum(np.zeros(len(self.v_n.x.array)), self.v_n.x.array)

        if self.param.mu == 1:

            with self.b_v.localForm() as loc_b_v:
                loc_b_v.set(0)
            fem.petsc.assemble_vector(self.b_v, self.linear_form_v)
            
            # Apply Dirichlet boundary condition to the vector
            fem.petsc.apply_lifting(self.b_v, [self.bilinear_form_v], [self.geometry.bc])
            self.b_v.ghostUpdate(addv=PETSc.InsertMode.ADD_VALUES, mode=PETSc.ScatterMode.REVERSE)
            fem.petsc.set_bc(self.b_v, self.geometry.bc)
        
            # Solve linear problem
            self.solver_v.solve(self.b_v, self.v_n.vector)
            self.v_n.x.scatter_forward()
        
            self.v_n.x.array[:] = np.maximum(np.zeros(len(self.v_n.x.array)), self.v_n.x.array)

    def plot(self, n: int) -> None:

        PATH = op.join(file_dir, 'Figures', self.param.example_name)
        if not op.exists(PATH):
            os.makedirs(PATH)
        
        t = n * self.param.dt
        file_name_plot_u = op.join(PATH, f'figure_u_{self.i}.xdmf')
        file_name_plot_v = op.join(PATH, f'figure_v_{self.i}.xdmf')

        if n == 0:
            with io.XDMFFile(self.geometry.domain.comm, file_name_plot_u, 'w') as xdmf:
                xdmf.write_mesh(self.geometry.domain)
                xdmf.write_function(self.u_n, t)

            with io.XDMFFile(self.geometry.domain.comm, file_name_plot_v, 'w') as xdmf:
                xdmf.write_mesh(self.geometry.domain)
                xdmf.write_function(self.v_n, t)

        if n > 0:
            with io.XDMFFile(self.geometry.domain.comm, file_name_plot_u, 'a') as xdmf:
                xdmf.write_function(self.u_n, t)

            with io.XDMFFile(self.geometry.domain.comm, file_name_plot_v, 'a') as xdmf:
                xdmf.write_function(self.v_n, t)

    def save_results(self) -> None:
        PATH = op.join(file_dir, 'Results', self.param.example_name)
        if not op.exists(PATH):
            os.makedirs(PATH)

        file_name_results = op.join(PATH, f'results_{self.i}_{self.j}_{self.k}.pickle')
        with open(file_name_results, 'wb') as file:
            pickle.dump(self.Counter_list, file)
            pickle.dump(self.Error_L2, file)
            pickle.dump(self.Error_grad, file)
            pickle.dump(self.Error_main, file)


def pme_M_scheme(solution: solution_class, parameters: general_param) -> None:

    # Ininitialize u_n and v_n.
    solution.initialize_sol()
    if parameters.save_plot:
        solution.plot(0)

    # Loop over time steps.
    total_steps = int(solution.param.final_time/solution.param.dt)
    for n in range(total_steps):

        # Solve for (u_n, w_n) using M-scheme.
        solution.update_uw()

        if solution.Counter_list[-1] == -1:
            print('Stopped example as did not converge')
            break

        # Solve for v_n.
        solution.update_v()

        if parameters.save_plot:
            solution.plot(n+1)

        print(f'Scheme is at {(n+1)/total_steps * 100}%. \r')
            
    if parameters.save_results:
        print('Saving results.')
        solution.save_results()
    