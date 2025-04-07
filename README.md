# bio_IVP_MCMC
Quantitative assessment of biological dynamics with aggregate data. Supporting code of the paper [*Quantitative assessment of biological dynamics with aggregate data*](https://arxiv.org/abs/2504.02581).

This implementation of the code does the following:

1. Takes in pregenerated synthetic data based on real *Prochlorococcus* experiment data and uses this information to initialize an MCMC sampler.
2. Runs tailor-made Hamiltonian Monte Carlo (HMC) and Multiplicative Elliptical Slice Sampling (MESS) samplers to generate samples from posterior distribution described in the corresponding paper.
3. Tracks samples and visualizes data fit as well as the log posterior values.
4. Saves size-controlled chain for further exploration.

# Instructions


## Provide data 

The **synthetic_data.mat** file shows the experiment information needed. The time series values (t) along with the max and min. The mean values of the samples at each time step (z), the variance at each time step (v), and the number of data points being aggregated (K). 

## Set hyper-parameter and constraint condtions
In  **chainer_init_params.m** the hyper-parameters of the model as well as setting for the sampler can be adjusted. Also the constraint type can be adjusted depending on the data you provide. 

## Set ODE model and solvers
**get_x.m** holds the ODE solver

## Run Code
Use **run_chainer.m** to set output directory and file name where the chain will be stored. This function ends with **run_expander** which is a script to make repeat chain updates

## Post-processing
Once the sampler has been run in the **/Runs** directory the resulting chain will be stored as a .mat file where it can loaded into whatever type of post-processing scripts you need. 