# Author: Robin Smeets 
# Email: robinsmeets99@gmail.com / r.k.h.smeets@uva.nl
# Institute: Korteweg-de Vries Institute for Mathematics - University of Amsterdam

'''
Python script for generating and saving simulations of the pme model.
'''
import argparse
import pme_utils as utils
import pme_config as config

def generate_plot(geometry_config: utils.general_param, 
                  solution_config: utils.solution_param, 
                  general_config: utils.general_param,
                  parsed_args: argparse.Namespace) -> None:

    geometry = utils.geometry_class(geometry_config)
    solution = utils.solution_class(geometry, solution_config)

    utils.pme_M_scheme(solution, general_config, parsed_args, print_progress= True)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-mv', '--movewsl', nargs='?', type=bool, const=True, default=False, 
                        help= 'If true, copies the .bp files to a specified folder (general_param.windows_dir) in Windows from WSL2')
    parsed_args, unknown = parser.parse_known_args()

    # # 1D test simulation
    # geometry_config = config.geometry_param_1D_simulation_pme
    # solution_config = config.solution_param_1D_simulation_pme
    # general_config = config.general_param_1D_simulation_pme
    # generate_plot(geometry_config, solution_config, general_config, parsed_args)

    # # 2D test simulation
    # geometry_config = config.geometry_param_2D_simulation_pme
    # solution_config = config.solution_param_2D_simulation_pme
    # general_config = config.general_param_2D_simulation_pme
    # generate_plot(geometry_config, solution_config, general_config, parsed_args)

