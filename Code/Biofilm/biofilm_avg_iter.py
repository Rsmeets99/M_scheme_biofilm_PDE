from inspect import getsourcefile
import os
from os import path as op
import pickle

import matplotlib.pyplot as plt
import numpy as np

import biofilm_utils as utils
import biofilm_config as config

file_name = op.abspath(getsourcefile(lambda:0))
file_dir = op.dirname(file_name)
 
def generate_data(dt_list: list[float], M_list: list[float], h_list: list[float]) -> None:
    geometry_config = config.geometry_param_1D_avg_iter_biofilm
    solution_config = config.solution_param_1D_avg_iter_biofilm
    general_config = config.general_param_avg_iter_biofilm

    for i, dt in enumerate(dt_list):
        for j, M_par in enumerate(M_list):
            for k, h in enumerate(h_list):
                geometry_config.h = h
                solution_config.dt = dt
                solution_config.M_par = M_par

                geometry = utils.geometry_class(geometry_config)
                solution = utils.solution_class(geometry, solution_config)

                solution.i, solution.j, solution.k = i, j, k

                utils.biofilm_M_scheme(solution, general_config)

def combine_data(dt_list: list[float], M_list: list[float], h_list: list[float]) -> None:
    solution_config = config.solution_param_1D_avg_iter_biofilm

    results_dir = op.join(file_dir, 'Results', 'Combined')
    if not op.exists(results_dir):
        os.makedirs(results_dir)
    name = solution_config.example_name
    example_dir = op.join(results_dir, name)
    if not op.exists(example_dir):
        os.makedirs(example_dir)

    for i, dt in enumerate(dt_list):
        path_name = op.join(example_dir, f'avg_iter_results_dt_{dt}.txt')
        with open(path_name, 'w', newline='') as file:
            file.write('h\t' + '\t'.join(f'{M}' for M in M_list) + '\n')

        temp_array = np.zeros((len(h_list), len(M_list)))
        
        for k, h in enumerate(h_list):
            for j, M_par in enumerate(M_list):

                file_name_results = op.join(example_dir, f'results_{i}_{j}_{k}.pickle')
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
                    

if __name__ == '__main__':

    dt_list = [10**(-1*(1+i/2)) for i in range(0,4)]
    h_list = [1/i for i in range(10,210,10)]
    hinv_list = [i for i in range(10,210,10)]
    M_list = [1e-7, 5e-5, 5e-4, 5e-3]

    generate_data(dt_list, M_list, h_list)
    combine_data(dt_list, M_list, h_list)




