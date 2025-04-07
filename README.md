# bio_IVP_MCMC
Quantitative assessment of biological dynamics with aggregate data. Supporting code of the paper [*Quantitative assessment of biological dynamics with aggregate data*](https://arxiv.org/abs/2504.02581).

This implementation of the code does the following:

1. Takes in pregenerated synthetic data based on real *Prochlorococcus* experiment data and uses this information to initialize an MCMC sampler.
2. Runs tailor made Hamiltonian Monte Carlo (HMC) and Multiplicative Elliptical Slice Sampling (MESS) samplers to generate samples from posterior distribution described in the corresponding paper.
3. Tracks samples and visualizes data fit as well as the log posterior values.
4. Saves size controlled chain for further exploration.

# Instructions


## Provide data 

## Set hyper-parameter and constraint condtions

## Set ODE model and solvers

## Run Code

## Post-processing