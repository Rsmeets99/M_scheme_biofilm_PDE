# General explanation of this repository
This repository contains all the code used to replicate the results within our paper [Robust time-discretisation and linearisation schemes for singular and degenerate evolution systems modelling biofilm growth](https://arxiv.org/abs/2404.00391). The results within the paper were generated with an older version of the code written during my master's thesis in FEniCSx v0.5.1, but for this paper the code has been completely rewritten to make it better readable, easier to customize and updated to the most recent stable version of FEniCSx (v0.7.3) as of paper submission. 

If one is interested in the older version for FEniCSx v0.5.1 or if one has any other questions or comments, please feel free to send me an e-mail (robinsmeets99@gmail.com).

# How to install
The easiest way to install the prerequisite packages (in my experience) is to first install miniconda or anaconda, and use that to create an environment that contains the packages. See also the official FEniCSx download site [link](https://fenicsproject.org/download/) (or for more details [link](https://github.com/FEniCS/dolfinx#installation))for a how-to. Furthermore, one requires numpy and matplotlib to generate the figures. This means that after installing miniconda or anaconda, one can use the commands

```
conda create -n fenicsx-env
conda activate fenicsx-env
conda install -c conda-forge fenics-dolfinx mpich pyvista numpy matploblib
```

to create an environment in which the code can be run. Note that if ones uses windows, the installation must be done through the Windows Subsystem for Linux (WSL/WSL2) and install Ubuntu. Alternatively, FEniCSx can be ran through Docker.

# How to run
Explanation on the different .py scripts is given within the README.md files of their respective directories.

# How to view simulations
The simulation data is stored with VTX in a .bp folder. Note that these can be quite big (on the order of gb for larger simulations), which is the reason I could not upload any premade simulations to GitHub (only accepts files smaller than 100 mb). These .bp files can be viewed within [ParaView](https://www.paraview.org/download/). You want to use the ADIOSVTX2READER to open the data. For 1D one needs the filter `plot data`, while for 2D one needs the filter `scale by scalar` to get a 3D plot. For 2D it might also be necessary to go to the properties of the solution and select coloring and then choose u_n instead of solid color to get the plot. Afterwards, one can then use the scale by scalar filter to make a 3D plot.

# How to cite
When using the results or code within this repository, we ask you kindly to cite our paper ([link to Arxiv](https://arxiv.org/abs/2404.00391)).
```
@article{SMSP24,
  author={Smeets, Robin and Mitra, Koondanibha and Pop, Sorin and Sonner, Stefanie},
  journal={arXiv preprint arXiv:2404.00391},
  title={Robust time-discretisation and linearisation schemes for singular and degenerate evolution systems modelling biofilm growth},
  year={2024}
}
```