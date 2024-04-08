This repository is a work in progress. Currently in the progress of completely refactoring the old code and making it compatible with the newest stable version of fenicsx. If one wants earlier acces to the code, please contact me for the old legacy code. 

# TODO

- [ ] Add plotting functions for `biofilm_avg_iter` and `biofilm_conv_rate`
  - [ ] Include fitting of straight line in log-log scale for `biofilm_conv_rate`
- [ ] Bug testing this new version of code
- [ ] Repeat the work for the porous medium case
  - [ ] Do not forget to add convergence to exact solution as dt to 0
- [ ] Add comments and docstrings to code
- [ ] Finish the `README.md` files

# General explanation of this repository

# How to install
Install miniconda/anaconda, create a new environment with fenicsx. Need to expand on this.

# How to run

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