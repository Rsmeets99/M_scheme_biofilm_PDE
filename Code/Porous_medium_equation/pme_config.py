
from pme_utils import geometry_param, solution_param, general_param

## Configurations for the average amount of iterations against mesh-size.

geometry_param_1D_avg_iter_pme_dict = {
    'dim' : 1,
    'h' : 0.1,
    'x_min' : -2,
    'x_max' : 2
}

solution_param_1D_avg_iter_pme_dict = {
    'example_name' : 'pme_avg_iter_1D',

    'start_time' : 0.3,
    'final_time' : 1.2,
    'dt' : 0.1,
    'gamma': 1/4,
    'M_par': 1e-7,

    'stop_crit': 1e-5,
    'allowed_iter' : 400,
    'full_iter' : False,

    'r' : 0.2,
    'height' : 0.9,
    'xlist' : [-0.3, 0.3],
    'ylist' : [-0.3, 0.3],

    'un_max_estimate': 0.993,
}

general_param_avg_iter_pme_dict = {
    'save_plot' : False,
    'save_results' : True
}

geometry_param_1D_avg_iter_pme = geometry_param(**geometry_param_1D_avg_iter_pme_dict)
solution_param_1D_avg_iter_pme = solution_param(**solution_param_1D_avg_iter_pme_dict)
general_param_avg_iter_pme = general_param(**general_param_avg_iter_pme_dict)

## Configurations for the convergence rate alpha against dt.

geometry_param_1D_conv_alpha_pme_dict = {
    'dim' : 1,
    'h' : 1e-4,
    'x_min' : -2,
    'x_max' : 2
}

solution_param_1D_conv_alpha_pme_dict = {
    'example_name' : 'pme_conv_alpha_1D',

    'start_time' : 0.3,
    'final_time' : 0.1,
    'dt' : 0.1,
    'gamma': 1/4,
    'M_par': 5e-4,

    'stop_crit': 1e-5,
    'allowed_iter' : 4,
    'full_iter' : True,

    'r' : 0.2,
    'height' : 0.9,
    'xlist' : [-0.3, 0.3],
    'ylist' : [-0.3, 0.3],

    'un_max_estimate': 0.993,
}

general_param_conv_alpha_pme_dict = {
    'save_plot' : False,
    'save_results' : True
}

geometry_param_1D_conv_alpha_pme = geometry_param(**geometry_param_1D_conv_alpha_pme_dict)
solution_param_1D_conv_alpha_pme = solution_param(**solution_param_1D_conv_alpha_pme_dict)
general_param_conv_alpha_pme = general_param(**general_param_conv_alpha_pme_dict)



if __name__ == '__main__':

    print(geometry_param_1D_avg_iter_pme)
    print(solution_param_1D_avg_iter_pme)
    print(general_param_avg_iter_pme)

    print(geometry_param_1D_conv_alpha_pme)
    print(solution_param_1D_conv_alpha_pme) 
    print(general_param_conv_alpha_pme) 
