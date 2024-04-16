# Author: Robin Smeets 
# Email: robinsmeets99@gmail.com / r.k.h.smeets@uva.nl
# Institute: Korteweg-de Vries Institute for Mathematics - University of Amsterdam

'''
Python script containing the configurations of all the different examples.
'''

from pme_utils import geometry_param, solution_param, general_param, experiment_param

## Configurations for the average amount of iterations against mesh-size.
geometry_param_1D_avg_iter_pme_dict = {
    'dim' : 1,
    'h' : 0.1,
    'x_min' : -2,
    'x_max' : 2
}

solution_param_1D_avg_iter_pme_dict = {
    'example_name' : 'pme_avg_iter_1D',      

    'start_time' : 0.5,       
    'final_time' : 1.1,     
    'dt' : 0.1,             
    'gamma' : 1/3,            
    'M_par' : 1e-3,            

    'm' : 4,               
    'beta' : 1,             
    'C' : 1/32,                

    'stop_crit' : 1e-5,        
    'allowed_iter' : 400,       
    'full_iter' : False,         
}

general_param_1D_avg_iter_pme_dict = {
    'save_plot' : False,
    'save_results' : True,
    'save_time_error' : False
}

experiment_param_1D_avg_iter_pme_dict = {
    'dt_list' : [10**(-1*(1+i/2)) for i in range(0,4)],
    'h_list' : [1/i for i in range(10,210,10)],
    'M_list' : [1e-7, 5e-5, 5e-4, 5e-3],
}

geometry_param_1D_avg_iter_pme = geometry_param(**geometry_param_1D_avg_iter_pme_dict)
solution_param_1D_avg_iter_pme = solution_param(**solution_param_1D_avg_iter_pme_dict)
general_param_1D_avg_iter_pme = general_param(**general_param_1D_avg_iter_pme_dict)
experiment_param_1D_avg_iter_pme = experiment_param(**experiment_param_1D_avg_iter_pme_dict)

## Configurations for the convergence rate alpha against dt.
geometry_param_1D_conv_alpha_pme_dict = {
    'dim' : 1,
    'h' : 1e-4,
    'x_min' : -2,
    'x_max' : 2
}

solution_param_1D_conv_alpha_pme_dict = {
    'example_name' : 'pme_conv_alpha_1D',

    'start_time' : 0.5,       
    'final_time' : 1.1,     
    'dt' : 0.1,             
    'gamma' : 1/3,            
    'M_par' : 1e-3,            

    'm' : 4,               
    'beta' : 1,             
    'C' : 1/32,                

    'stop_crit' : 1e-5,        
    'allowed_iter' : 4,       
    'full_iter' : True,
}

general_param_1D_conv_alpha_pme_dict = {
    'save_plot' : False,
    'save_results' : True,
    'save_time_error' : False
}

experiment_param_1D_conv_alpha_pme_dict = {
    'dt_list' : [10**(-i/4) for i in range(4,11)],
}

geometry_param_1D_conv_alpha_pme = geometry_param(**geometry_param_1D_conv_alpha_pme_dict)
solution_param_1D_conv_alpha_pme = solution_param(**solution_param_1D_conv_alpha_pme_dict)
general_param_1D_conv_alpha_pme = general_param(**general_param_1D_conv_alpha_pme_dict)
experiment_param_1D_conv_alpha_pme = experiment_param(**experiment_param_1D_conv_alpha_pme_dict)

## Configurations for the dt convergence rate to the exact solution.
geometry_param_1D_conv_time_pme_dict = {
    'dim' : 1,
    'h' : 1e-4,
    'x_min' : -2,
    'x_max' : 2
}

solution_param_1D_conv_time_pme_dict = {
    'example_name' : 'pme_conv_time_1D',

    'start_time' : 0.5,       
    'final_time' : 1,     
    'dt' : 0.1,             
    'gamma' : 1/3,            
    'M_par' : 1e-3,            

    'm' : 4,               
    'beta' : 1,             
    'C' : 1/32,                

    'stop_crit' : 1e-7,        
    'allowed_iter' : 400,       
    'full_iter' : False,
}

general_param_1D_conv_time_pme_dict = {
    'save_plot' : False,
    'save_results' : True,
    'save_time_error' : True
}

experiment_param_1D_conv_time_pme_dict = {
    'dt_list' : [10**(-i/4) for i in range(4,13)],
}

geometry_param_1D_conv_time_pme = geometry_param(**geometry_param_1D_conv_time_pme_dict)
solution_param_1D_conv_time_pme = solution_param(**solution_param_1D_conv_time_pme_dict)
general_param_1D_conv_time_pme = general_param(**general_param_1D_conv_time_pme_dict)
experiment_param_1D_conv_time_pme = experiment_param(**experiment_param_1D_conv_time_pme_dict)

## Configurations for 1D test simulation
geometry_param_1D_simulation_pme_dict = {
    'dim' : 1,
    'h' : 1e-4,
    'x_min' : -2,
    'x_max' : 2
}

solution_param_1D_simulation_pme_dict = {
    'example_name' : 'pme_simulation_1D',

    'start_time' : 0.5,       
    'final_time' : 1.1,     
    'dt' : 0.01,             
    'gamma' : 1/3,            
    'M_par' : 1e-3,            

    'm' : 4,               
    'beta' : 1,             
    'C' : 1/32,                

    'stop_crit' : 1e-5,        
    'allowed_iter' : 400,       
    'full_iter' : False,
}

general_param_1D_simulation_pme_dict = {
    'save_plot' : True,
    'save_results' : True,
    'save_time_error' : False,
    'windows_dir' : '/mnt/c/Users/rsmeets/OneDrive - UvA/Documenten/Master Thesis/Simulation bp files'
}

geometry_param_1D_simulation_pme = geometry_param(**geometry_param_1D_simulation_pme_dict)
solution_param_1D_simulation_pme = solution_param(**solution_param_1D_simulation_pme_dict)
general_param_1D_simulation_pme = general_param(**general_param_1D_simulation_pme_dict)

## Configurations for 2D simulations of pme PDE-PDE.
geometry_param_2D_simulation_pme_dict = {
    'dim' : 2,
    'h' : 0.02,
    'x_min' : -2,
    'x_max' : 2
}

solution_param_2D_simulation_pme_dict = {
    'example_name' : 'pme_simulation_2D',

    'start_time' : 0.5,       
    'final_time' : 1.1,     
    'dt' : 0.01,             
    'gamma' : 1/3,            
    'M_par' : 1e-3,            

    'm' : 4,               
    'beta' : 1,             
    'C' : 1/32,                

    'stop_crit' : 1e-5,        
    'allowed_iter' : 400,       
    'full_iter' : False,
}

general_param_2D_simulation_pme_dict = {
    'save_plot' : True,
    'save_results' : True,
    'save_time_error' : False,
    'windows_dir' : '/mnt/c/Users/rsmeets/OneDrive - UvA/Documenten/Master Thesis/Simulation bp files'
}

geometry_param_2D_simulation_pme = geometry_param(**geometry_param_2D_simulation_pme_dict)
solution_param_2D_simulation_pme = solution_param(**solution_param_2D_simulation_pme_dict)
general_param_2D_simulation_pme = general_param(**general_param_2D_simulation_pme_dict)


if __name__ == '__main__':

    print(geometry_param_1D_avg_iter_pme)
    print(solution_param_1D_avg_iter_pme)
    print(general_param_1D_avg_iter_pme)
    print(experiment_param_1D_avg_iter_pme)

    print(geometry_param_1D_conv_alpha_pme)
    print(solution_param_1D_conv_alpha_pme) 
    print(general_param_1D_conv_alpha_pme)
    print(experiment_param_1D_conv_alpha_pme) 

    print(geometry_param_1D_conv_time_pme)
    print(solution_param_1D_conv_time_pme)
    print(general_param_1D_conv_time_pme)
    print(experiment_param_1D_conv_time_pme)

    print(geometry_param_1D_simulation_pme) 
    print(solution_param_1D_simulation_pme) 
    print(general_param_1D_simulation_pme) 

    print(geometry_param_2D_simulation_pme) 
    print(solution_param_2D_simulation_pme)
    print(general_param_2D_simulation_pme) 
