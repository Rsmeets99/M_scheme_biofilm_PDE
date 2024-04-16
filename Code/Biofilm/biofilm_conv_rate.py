# Author: Robin Smeets 
# Email: robinsmeets99@gmail.com / r.k.h.smeets@uva.nl
# Institute: Korteweg-de Vries Institute for Mathematics - University of Amsterdam

'''
Python script for computing the convergence rate against the time-step size.
'''

import argparse
import datetime
from inspect import getsourcefile
import os
from os import path as op
import pickle
import time

import numpy as np
import numpy.typing as npt
import matplotlib.pyplot as plt

import biofilm_utils as utils
import biofilm_config as config

file_name = op.abspath(getsourcefile(lambda:0))
file_dir = op.dirname(file_name)

def generate_data(dt_list: list[float], 
                  geometry_config: utils.geometry_param, 
                  solution_config: utils.solution_param, 
                  general_config: utils.general_param,
                  parsed_args: argparse.Namespace) -> None:
    '''Generate the data for plotting the convergence rate against time-step size.'''

    geometry = utils.geometry_class(geometry_config)

    time_start_readable_all = datetime.datetime.fromtimestamp(time.time())
    print(f'Start generating the data for example {solution_config.example_name}. Current time: {time_start_readable_all} \n')

    for i, dt in enumerate(dt_list):
        print(f'\rGenerating data for dt = {dt:.2E}', end='')
        solution_config.dt = dt
        solution_config.final_time = dt
        solution = utils.solution_class(geometry, solution_config)
        solution.i = i

        utils.biofilm_M_scheme(solution, general_config, parsed_args)

    time_end_readable_all = datetime.datetime.fromtimestamp(time.time())
    print(f'\nFinished generating the data for example {solution_config.example_name}. Current time: {time_end_readable_all}')
    print(f'It took {time_end_readable_all-time_start_readable_all}')

def combine_data(dt_list: list[float], 
                 solution_config: utils.solution_param) -> None:
    '''Generate a .txt file with the relevant generated data.'''

    name = solution_config.example_name
        
    results_dir = op.join(file_dir, 'Results')
    if not op.exists(results_dir):
        os.makedirs(results_dir)
    example_dir = op.join(file_dir, 'Processed', name)
    if not op.exists(example_dir):
        os.makedirs(example_dir)

    path_name = op.join(example_dir, f'conv_rate_results.txt')
    with open(path_name, 'w', newline='') as file:
        file.write('dt\trate\n')

    for i, dt in enumerate(dt_list):

        file_name_results = op.join(results_dir, name, f'results_{i}_0_0.pickle')
        data = []
        with open(file_name_results, 'rb') as openfile:
            while True:
                try:
                    data.append(pickle.load(openfile))
                except EOFError:
                    break
        
        Error = data[3][0]
        geometric_mean_fraction = (Error[-1]/Error[0])**(1/len(Error))
        with open(path_name, 'a', newline='') as file:
            file.write(f'{dt}\t{geometric_mean_fraction}\n')

def linear_log_fit(dt_data: npt.ArrayLike, conv_data: npt.ArrayLike, cutoff : int = 0) -> tuple[npt.ArrayLike, float]:
    '''Fit a linear line through the log of the data and output the corresponding exponential curve'''

    log_dt_data = np.log10(dt_data)
    log_conv_data = np.log10(conv_data)

    a, b = np.polyfit(log_dt_data[cutoff:], log_conv_data[cutoff:], 1)
    line_data = 10**(b) * dt_data**(a)

    return line_data, a


def plot_data(geometry_config: utils.geometry_param, 
              solution_config: utils.solution_param) -> None:
    '''Plot the figures from the combined data. '''
    
    name = solution_config.example_name

    results_dir = op.join(file_dir, 'Processed')
    example_dir = op.join(results_dir, name)

    path_name = op.join(example_dir, f'conv_rate_results.txt')
    path_name_fig = op.join(example_dir, f'conv_rate_results_figure.png')

    data = []
    with open(path_name) as file:
        for line in file:
            row = line.split('\t')
            data.append(row)

    true_data = np.asarray(data[1:], dtype=float)

    dt_data = true_data[:,0]
    conv_data = true_data[:,1]

    line_fit, a = linear_log_fit(dt_data, conv_data)

    plt.plot(np.log10(dt_data), np.log10(conv_data), 'r.', label = '$\\alpha$')
    plt.plot(np.log10(dt_data), np.log10(line_fit), 'k--', label = f'slope = {a:.2f}')

    plt.xlabel('Log10 time step size $\\tau$')
    plt.ylabel('Log10 convergence rate $\\alpha$')
    plt.title(f'Convergence rate, $h = {geometry_config.h:.2E}$, $\\gamma = {solution_config.gamma:.2f}$, {geometry_config.dim:.0f}D, $\\mu = {solution_config.mu:.0f}$')
    plt.grid()
    plt.legend()
    plt.savefig(path_name_fig)
    plt.clf()

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-mv', '--movewsl', nargs='?', type=bool, const=True, default=False, 
                        help= 'If true, copies the .bp files to a specified folder (general_param.windows_dir) in Windows from WSL2.')
    parser.add_argument('-gd', '--generatedata', nargs='?', type=bool, const=True, default=False, 
                        help= 'If true, generates new data, overwriting existing data.')
    parser.add_argument('-cd', '--combinedata', nargs='?', type=bool, const=True, default=False, 
                        help= 'If true, combines new data, overwriting existing combined data.')
    parser.add_argument('-pd', '--plotdata', nargs='?', type=bool, const=False, default=True, 
                        help= 'If true, plots combined data, overwriting existing plots.')
    parsed_args, unknown = parser.parse_known_args()

    dt_list = config.experiment_param_1D_conv_alpha_biofilm.dt_list

    geometry_config = config.geometry_param_1D_conv_alpha_biofilm
    solution_config = config.solution_param_1D_conv_alpha_biofilm
    general_config = config.general_param_1D_conv_alpha_biofilm

    if parsed_args.generatedata:
        generate_data(dt_list, geometry_config, solution_config, general_config, parsed_args)
    if parsed_args.combinedata:
        combine_data(dt_list, solution_config)
    if parsed_args.plotdata:
        plot_data(geometry_config, solution_config)
