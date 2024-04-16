# Author: Robin Smeets 
# Email: robinsmeets99@gmail.com / r.k.h.smeets@uva.nl
# Institute: Korteweg-de Vries Institute for Mathematics - University of Amsterdam

'''
Python script for computing the average amount of required iterations for different time steps against the mesh size.
'''

import argparse
import datetime
from inspect import getsourcefile
import os
from os import path as op
import pickle
import time

import matplotlib.pyplot as plt
import numpy as np

import biofilm_utils as utils
import biofilm_config as config

file_name = op.abspath(getsourcefile(lambda:0))
file_dir = op.dirname(file_name)
 
def generate_data(dt_list: list[float], 
                  M_list: list[float], 
                  h_list: list[float], 
                  geometry_config: utils.geometry_param, 
                  solution_config: utils.solution_param, 
                  general_config: utils.general_param,
                  parsed_args: argparse.Namespace) -> None:
    '''Generate the data for plotting the average amount of iterations.'''

    time_start_readable_all = datetime.datetime.fromtimestamp(time.time())
    print(f'Start generating the data for example {solution_config.example_name}.\nStarting time: {time_start_readable_all}\n')

    for i, dt in enumerate(dt_list):
        time_start = datetime.datetime.fromtimestamp(time.time())
        print(f'Start generating data for dt = {dt:.2E}, time is {time_start}')
        for k, h in enumerate(h_list):
            print(f'\rCurrent h is {h:.2E}', end='')
            geometry_config.h = h
            geometry = utils.geometry_class(geometry_config)

            for j, M_par in enumerate(M_list):
                solution_config.dt = dt
                solution_config.M_par = M_par

                solution = utils.solution_class(geometry, solution_config)
                solution.i, solution.j, solution.k = i, j, k

                utils.biofilm_M_scheme(solution, general_config, parsed_args)
        
        time_end = datetime.datetime.fromtimestamp(time.time())
        print(f'\rIt took {time_end - time_start} to get the results for dt = {dt:.2E}\n')
    
    time_end_readable_all = datetime.datetime.fromtimestamp(time.time())
    print(f'Finished generating the data for example {solution_config.example_name}.\nFinishing time: {time_end_readable_all}')
    print(f'It took {time_end_readable_all-time_start_readable_all}')

def combine_data(dt_list: list[float], 
                 M_list: list[float], 
                 h_list: list[float],
                 solution_config: utils.solution_param) -> None:
    '''Generate a .txt file with the combined relevant generated data.'''

    name = solution_config.example_name

    results_dir = op.join(file_dir, 'Results')
    if not op.exists(results_dir):
        os.makedirs(results_dir)
    example_dir = op.join(file_dir, 'Processed', name)
    if not op.exists(example_dir):
        os.makedirs(example_dir)

    for i, dt in enumerate(dt_list):
        path_name = op.join(example_dir, f'avg_iter_results_dt_{dt:.2E}.txt')
        with open(path_name, 'w', newline='') as file:
            file.write('h\t' + '\t'.join(f'{M}' for M in M_list) + '\n')

        temp_array = np.zeros((len(h_list), len(M_list)))
        
        for k, h in enumerate(h_list):
            for j, M_par in enumerate(M_list):

                file_name_results = op.join(results_dir, name, f'results_{i}_{j}_{k}.pickle')
                data = []
                with open(file_name_results, 'rb') as openfile:
                    while True:
                        try:
                            data.append(pickle.load(openfile))
                        except EOFError:
                            break
                
                iterations = data[0]
                if iterations[-1] == -1:
                    temp_array[k,j] = None
                else:
                    avg_iterations = np.mean(iterations)
                    temp_array[k,j] = avg_iterations
            
            with open(path_name, 'a', newline='') as file:
                file.write(f'{h}\t' + '\t'.join(f'{avg}' for avg in temp_array[k,:]) + '\n')


def plot_data(dt_list: list[float], 
              solution_config: utils.solution_param) -> None:
    '''Plot the figures from the combined data.'''

    name = solution_config.example_name

    results_dir = op.join(file_dir, 'Processed')
    example_dir = op.join(results_dir, name)

    for dt in dt_list:
        
        path_name = op.join(example_dir, f'avg_iter_results_dt_{dt:.2E}.txt')
        path_name_fig = op.join(example_dir, f'avg_iter_results_dt_{dt:.2E}_figure.png')

        data = []
        with open(path_name) as file:
            for line in file:
                row = line.split('\t')
                data.append(row)

        labels = data[0]
        true_data = np.asarray(data[1:], dtype=float)

        for i in range(len(labels)-1):
            plt.plot(1/true_data[:,0], true_data[:,i+1], '--', marker = '.', markersize=10, label = f'M_par = {float(labels[i+1]):.2E}')

        plt.xlabel('Mesh size 1/h')
        plt.ylabel('Iterations')
        plt.title(f'Average amount of iterations, $\\tau = {dt:.1E}$ and $T = 1.2$')
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

    # dt_list = config.experiment_param_1D_avg_iter_biofilm.dt_list
    # h_list = config.experiment_param_1D_avg_iter_biofilm.h_list
    # M_list = config.experiment_param_1D_avg_iter_biofilm.M_list

    # Smaller parameters for debugging
    dt_list = [10**(-1*(1+i/2)) for i in range(0,2)]
    h_list = [1/i for i in range(10,60,10)]
    M_list = [1e-7, 5e-4]

    geometry_config = config.geometry_param_1D_avg_iter_biofilm
    solution_config = config.solution_param_1D_avg_iter_biofilm
    general_config = config.general_param_1D_avg_iter_biofilm

    if parsed_args.generatedata:
        generate_data(dt_list, M_list, h_list, geometry_config, solution_config, general_config, parsed_args)
    if parsed_args.combinedata:
        combine_data(dt_list, M_list, h_list, solution_config)
    if parsed_args.plotdata:
        plot_data(dt_list, solution_config)
