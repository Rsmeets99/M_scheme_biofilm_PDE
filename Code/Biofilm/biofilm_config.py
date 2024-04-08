
from dataclasses import dataclass

from biofilm_utils import geometry_param, solution_param, general_param

## Configurations for the average amount of iterations against mesh-size.

geometry_param_1D_avg_iter_biofilm_dict = {
    'dim' : 1,
    'h' : 0.1,
    'x_min' : -2,
    'x_max' : 2
}

solution_param_1D_avg_iter_biofilm_dict = {
    'example_name' : 'biofilm_avg_iter_1D',

    'final_time' : 1.2,
    'dt' : 0.1,
    'gamma': 1/4,
    'M_par': 1e-7,
    'mu': 0,

    'stop_crit': 1e-5,
    'allowed_iter' : 400,
    'full_iter' : False,

    'r' : 0.2,
    'height' : 0.9,
    'xlist' : [-0.3, 0.3],
    'ylist' : [-0.3, 0.3],

    'un_max_estimate': 0.993,

    'kappa_1' : 0.4,
    'kappa_2' : 0.01,
    'kappa_3' : 1,
    'kappa_4' : 0.42,
    'delta_1' : 1*1e-6,
    'delta_2' : 0.2
}

general_param_avg_iter_biofilm_dict = {
    'save_plot' : False,
    'save_results' : True
}

geometry_param_1D_avg_iter_biofilm = geometry_param(**geometry_param_1D_avg_iter_biofilm_dict)
solution_param_1D_avg_iter_biofilm = solution_param(**solution_param_1D_avg_iter_biofilm_dict)
general_param_avg_iter_biofilm = general_param(**general_param_avg_iter_biofilm_dict)

## Configurations for the convergence rate alpha against dt.

geometry_param_1D_conv_alpha_biofilm_dict = {
    'dim' : 1,
    'h' : 1e-4,
    'x_min' : -2,
    'x_max' : 2
}

solution_param_1D_conv_alpha_biofilm_dict = {
    'example_name' : 'biofilm_conv_alpha_1D',

    'final_time' : 0.1,
    'dt' : 0.1,
    'gamma': 1/4,
    'M_par': 5e-4,
    'mu': 0,

    'stop_crit': 1e-5,
    'allowed_iter' : 3,
    'full_iter' : True,

    'r' : 0.2,
    'height' : 0.9,
    'xlist' : [-0.3, 0.3],
    'ylist' : [-0.3, 0.3],

    'un_max_estimate': 0.993,

    'kappa_1' : 0.4,
    'kappa_2' : 0.01,
    'kappa_3' : 1,
    'kappa_4' : 0.42,
    'delta_1' : 1*1e-6,
    'delta_2' : 0.2
}

general_param_conv_alpha_biofilm_dict = {
    'save_plot' : False,
    'save_results' : True
}

geometry_param_1D_conv_alpha_biofilm = geometry_param(**geometry_param_1D_conv_alpha_biofilm_dict)
solution_param_1D_conv_alpha_biofilm = solution_param(**solution_param_1D_conv_alpha_biofilm_dict)
general_param_conv_alpha_biofilm = general_param(**general_param_conv_alpha_biofilm_dict)

## Configurations for 2D simulations of biofilm PDE-PDE.

geometry_param_2D_PDEPDE_biofilm_dict = {
    'dim' : 2,
    'h' : 0.02,
    'x_min' : -2,
    'x_max' : 2
}

solution_param_2D_PDEPDE_biofilm_dict = {
    'example_name' : 'biofilm_plot_PDEPDE_2D',

    'final_time' : 0.1,
    'dt' : 0.01,
    'gamma': 1/4,
    'M_par': 1e-2,
    'mu': 1,

    'stop_crit': 1e-5,
    'allowed_iter' : 3,
    'full_iter' : True,

    'r' : 0.2,
    'height' : 0.9,
    'xlist' : [-0.3, 0.3],
    'ylist' : [-0.3, 0.3],

    'un_max_estimate': 0.99,

    'kappa_1' : 5,
    'kappa_2' : 0.01,
    'kappa_3' : 1,
    'kappa_4' : 0.42,
    'delta_1' : 5*1e-6,
    'delta_2' : 0.2
}

general_param_2D_PDEPDE_biofilm_dict = {
    'save_plot' : True,
    'save_results' : False
}

geometry_param_2D_PDEPDE_biofilm = geometry_param(**geometry_param_2D_PDEPDE_biofilm_dict)
solution_param_2D_PDEPDE_biofilm = solution_param(**solution_param_2D_PDEPDE_biofilm_dict)
general_param_2D_PDEPDE_biofilm = general_param(**general_param_conv_alpha_biofilm_dict)

## Configurations for 2D simulations of biofilm PDE-ODE.


geometry_param_2D_PDEODE_biofilm_dict = {
    'dim' : 2,
    'h' : 0.02,
    'x_min' : -2,
    'x_max' : 2
}

solution_param_2D_PDEODE_biofilm_dict = {
    'example_name' : 'biofilm_plot_PDEODE_2D',

    'final_time' : 0.1,
    'dt' : 0.01,
    'gamma': 1/4,
    'M_par': 1e-2,
    'mu': 0,

    'stop_crit': 1e-5,
    'allowed_iter' : 3,
    'full_iter' : True,

    'r' : 0.2,
    'height' : 0.9,
    'xlist' : [-0.3, 0.3],
    'ylist' : [-0.3, 0.3],

    'un_max_estimate': 0.99,

    'kappa_1' : 0.8,
    'kappa_2' : 0.01,
    'kappa_3' : 1,
    'kappa_4' : 0.42,
    'delta_1' : 8*1e-6,
    'delta_2' : 0.2
}

general_param_2D_PDEODE_biofilm_dict = {
    'save_plot' : True,
    'save_results' : False
}

geometry_param_2D_PDEODE_biofilm = geometry_param(**geometry_param_2D_PDEODE_biofilm_dict)
solution_param_2D_PDEODE_biofilm = solution_param(**solution_param_2D_PDEODE_biofilm_dict)
general_param_2D_PDEODE_biofilm = general_param(**general_param_2D_PDEODE_biofilm_dict)


if __name__ == '__main__':

    print(geometry_param_1D_avg_iter_biofilm)
    print(solution_param_1D_avg_iter_biofilm)
    print(general_param_avg_iter_biofilm)

    print(geometry_param_1D_conv_alpha_biofilm)
    print(solution_param_1D_conv_alpha_biofilm) 
    print(general_param_conv_alpha_biofilm) 

    print(geometry_param_2D_PDEPDE_biofilm) 
    print(solution_param_2D_PDEPDE_biofilm)
    print(general_param_2D_PDEPDE_biofilm) 

    print(geometry_param_2D_PDEODE_biofilm) 
    print(solution_param_2D_PDEODE_biofilm) 
    print(general_param_2D_PDEODE_biofilm)
