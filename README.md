# bio_IVP_MCMC
Quantitative assessment of biological dynamics with aggregate data. Supporting code of the paper [*Quantitative assessment of biological dynamics with aggregate data*](https://doi.org/10.48550/arXiv.2504.02581).

This source code does the following:

1. Loads pregenerated synthetic data based on real *Prochlorococcus* laboratory data and initializes an MCMC chain.
2. Runs tailor-made Hamiltonian Monte Carlo (cHMC) and Multiplicative Elliptical Slice Sampling (mESS) samplers to generate MCMC samples from the posterior distribution.
3. Tracks MCMC samples and produces fit visualization.
4. Saves MCMC chain for further exploration.

# Instructions


## Provide data 

The **example_data.mat** contains data properly formated.

## Set hyper-parameter and constraint condtions
In  **chainer_init_params.m** the hyper-parameters of the model as well as setting for the sampler can be adjusted. Also the constraint type can be adjusted depending on the data you provide. 

## Set ODE model and solvers
**get_x.m** implements the IVP solver

## Run code
Use **run_chainer.m** to set file name where the chain will be stored, initialize and expand an MCMC chain, and obtain basic posterior results.
