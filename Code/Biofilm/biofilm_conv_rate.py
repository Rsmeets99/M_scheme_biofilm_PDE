from inspect import getsourcefile
import os
from os import path as op
import pickle

import matplotlib.pyplot as plt

import biofilm_utils as utils
import biofilm_config as config

file_name = op.abspath(getsourcefile(lambda:0))
file_dir = op.dirname(file_name)

def generate_data(dt_list: list[float]) -> None:
    geometry_config = config.geometry_param_1D_conv_alpha_biofilm
    solution_config = config.solution_param_1D_conv_alpha_biofilm
    general_config = config.general_param_conv_alpha_biofilm

    geometry = utils.geometry_class(geometry_config)

    for i, dt in enumerate(dt_list):
        
        solution_config.dt = dt
        solution = utils.solution_class(geometry, solution_config)
        solution.i = i

        utils.biofilm_M_scheme(solution, general_config)


def combine_data(dt_list: list[float]) -> None:
    solution_config = config.solution_param_1D_conv_alpha_biofilm
    
    results_dir = op.join(file_dir, 'Results', 'Combined')
    if not op.exists(results_dir):
        os.makedirs(results_dir)
    name = solution_config.example_name
    example_dir = op.join(results_dir, name)
    if not op.exists(example_dir):
        os.makedirs(example_dir)

    path_name = op.join(example_dir, f'conv_rate_results.txt')
    with open(path_name, 'w', newline='') as file:
        file.write('dt\trate\n')

    for i, dt in enumerate(dt_list):

        file_name_results = op.join(example_dir, f'results_{i}_0_0.pickle')
        data = []
        with open(file_name_results, 'rb') as openfile:
            while True:
                try:
                    data.append(pickle.load(openfile))
                except EOFError:
                    break

        Error = data[3]
        geometric_mean_fraction = (Error[-1]/Error[0])**(1/len(Error))
        with open(path_name, 'a', newline='') as file:
            file.write(f'{dt}\t{geometric_mean_fraction}\n')


if __name__ == '__main__':

    dt_list = [10**(-i/4) for i in range(4,11)]

    generate_data(dt_list)
    combine_data(dt_list)
