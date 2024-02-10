# Files for Numerical Models for Random Performance Functions/Bimatrix Games/Replicator Equation

## Bimatrix Games
* Bimatrix games are turned into performance functions by means of 'bimatrix_Moran_sampler.m' given the payout matrices in the paper.
* For PD and stag, we did a spline on the log-odds and turned them into performance functions. Those functions are 'performance_pd_moran_24_individuals.m' and 'performance_stag_moran_24_individuals.m'. Since chicken was not smooth, instead we interpolated the function using 'Moran_performance_function_interp.m' on the data file ('log_odds_moran_chicken_10000_midway.mat').
* Gradient and Hessians are in the Hessian_for_performance_[game]_moran files.
* Then, run 'Evolution_Test_moran_softmax.m' with those files in the folder. 'PlotResults_Moran' plots the results.

## Replicator Equation
* Everything is in 'Evolution_Test_Switching_Analytic_Rankings.m' file.

## Random Performance Functions
* The different performance functions we tested is 'example_performance_6.m', and its gradient and Hessian are in 'Hessian_for_example_performance_6.m'.
* Once there, run the file 'Evolution_Test_RPF_softmax.m'. The plot functions are in 'PlotResults_RPFSM.m'.
