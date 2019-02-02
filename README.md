# Efficient Bayesian shape-constrained function estimation

R codes and function files for implementing the simulations and real data analysis addressed in the paper.

Description of the associated R files:

1. all_functions.R file contains all the required sub functions to incorporate Algorithm 1 of the paper, including functions like samp.WC    (draws a sample based on sampling scheme developed in Wood and Chan [1994]), ESS (draws a sample using elliptical slice sampler by        Murray et al. [2010]), nu.MH2 (draws posterior samples on the hyperparameters using Metropolis–Hastings and inverse Cholesky factor).      One can find detailed description of each functions in the file.
2. inv_chol.cpp file executes computations on inverse Cholesky factor utilizing Durbin’s recursion.
3. ESS_functions.R consists of three main functions that implement "shape restricted function estimation" for fixed hyperparameters. Three    main functions are separated based on the nature of the function fit: mon.inc.ESS (for monotone increasing), mon.dec.ESS (for monotone    decreasing) and con.ESS (for convex) functions take response y and covariate x as inputs, form the basis matrix, run the Gibbs sampler    and returns posterior samples on the parameters and function estimates along with 95% CI.
4. ESS_hyp_functions.R consists of three main functions that include hyperparameter updates for monotone increasing (mon.inc.ESS.hyp),        monotone decreasing (mon.dec.ESS.hyp) and convex (con.ESS.hyp) functions. These functions differ from the above three only by the          hyperparameter update steps.
5. ESsamp_hyp.R gives the code for Figure-2 (Section-4 : Simulation Results) of the paper.
6. cars_data_code.R is related to Section-5 (Application on a real data set) results on application of Algorithm 1. This code calculates      MSPE of one split of the 10 splits. Details can be found in the paper.
7. cars-mbart.csv contains car price data considered in Section-5 of the paper.

For more details on the codes or the functions, refer to the associated R files.
