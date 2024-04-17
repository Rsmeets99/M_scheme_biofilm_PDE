# Set-up porous medium equation folder
This folder contains all the code for the porous medium equation (pme) simulations and results as seen in the paper.

- [Set-up porous medium equation folder](#set-up-porous-medium-equation-folder)
- [General comments on running the code](#general-comments-on-running-the-code)
- [Processed](#processed)
- [Results](#results)
- [pme\_avg\_iter.py](#pme_avg_iterpy)
- [pme\_config.py](#pme_configpy)
- [pme\_conv\_rate.py](#pme_conv_ratepy)
- [pme\_simulation.py](#pme_simulationpy)
- [pme\_time\_conv.py](#pme_time_convpy)
- [pme\_utils.py](#pme_utilspy)

# General comments on running the code
Each of the python scripts has multiple flags for running the code: -gd, -cd, -pd, -mv
The simulation script only has -mv.

```
parser.add_argument('-mv', '--movewsl', nargs='?', type=bool, const=True, default=False, 
                        help= 'If true, copies the .bp files to a specified folder (general_param.windows_dir) in Windows from WSL2.')
parser.add_argument('-gd', '--generatedata', nargs='?', type=bool, const=True, default=False, 
                    help= 'If true, generates new data, overwriting existing data.')
parser.add_argument('-cd', '--combinedata', nargs='?', type=bool, const=True, default=False, 
                    help= 'If true, combines new data, overwriting existing combined data.')
parser.add_argument('-pd', '--plotdata', nargs='?', type=bool, const=False, default=True, 
                    help= 'If true, plots combined data, overwriting existing plots.')
```

The reason for -gd, -cd and -pd, is that one does not accidentally overwrite existing data, while the -mv flag exists as I have ParaView installed on Windows, while I run my code through WSL2, so I needed to transfer the .bp files to Windows.

# Processed
This folder contains the processed results of each of the python scripts: a .png plot of the wanted data as well as the cleaned data in a .txt file.

# Results
This folder contains all the generated .pickle files in running each of the python scripts. The relevant data in these pickle files is processed and saved within the `Processed` folder.

# pme_avg_iter.py
This python script generates and plots the data for the average amount of required iterations for convergence over some time span for our pme.

# pme_config.py
This python script contains all the parameters for each of the experiments (e.g. mesh size h, time step size dt, final time, etc.)

# pme_conv_rate.py
This python script generates and plots the data of the convergence rate $\alpha$ of the M-scheme for our pme against different time step sizes dt.

# pme_simulation.py
This python script generates simulations of the pme and saves them as .bp files within the `Simulation` folder. This folder is currently not there, as the size of these .bp files can easily reach the orders of gb which meant I could not upload them to GitHub. But if one were to run this file, one will find them within the then newly created `Simulation` folder.

# pme_time_conv.py
This python script generates and plots the data of the error between the numerical solution and the exact solution of our pme against different time step sizes dt.

# pme_utils.py
This python script contains all the mathematics. Here the pme functions are defined, the dataclasses are defined, a geometry class is defined which contains all the mesh related functions and a solution class is defined which has the M-scheme defined within.