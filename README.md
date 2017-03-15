# Summary

This repository contains all the files necessary to reproduce the
paper *[Expected Values Estimated via Mean-Field Approximation are
1/N-Accurate](mf_rate_convergence.pdf)* that has been accepted for at
ACM SIGMETRICS conference 2017.

Full reference : Nicolas Gast. Expected Values Estimated via
Mean-Field Approximation are 1/N-Accurate. Proceedings of ACM
SIGMETRICS conference 2017. Champaign-Urbana, USA, June 6-8, 2017
(SIGMETRICS’17).


# How to compile

To compile the paper: 
```{sh}
make paper 
```

The directories **jsq2_simulate** and **simu** contains the files
necessary to generate the figures.
* jsq2_simulate contains:
  * a file "two_choice_simulator.cc" that is a C++ implementation of a simulator of the two-choice model
  * an "ipython"-notebook that is used to generate simulation data and the figures
* simu contains the files to generates the figures of the small 1 dimensional example of Section 4. 

To compile the simulator and re-generate all figures, just type 
```{sh}
make simulator
make simulations
```
Note that the directory `jsq2_simulate/results` and `jsq2_simulate/traj` contain simulation results that are not re-generated each time (in order to save time). 
