# [Sampling of graph signals via randomized local aggregations](https://ieeexplore.ieee.org/abstract/document/8457237/)

BibTex reference:

@ARTICLE{Valsesia2019Sampling,
author={Diego {Valsesia} and Giulia {Fracastoro} and Enrico {Magli}},
journal={IEEE Transactions on Signal and Information Processing over Networks},
title={Sampling of Graph Signals via Randomized Local Aggregations},
year={2019},
volume={5},
number={2},
pages={348-359},
keywords={Signal processing;Compressed sensing;Stability analysis;Transforms;Information processing;Laplace equations;Sparse matrices;Graph signal processing;sampling;random projections;compressed sensing},
doi={10.1109/TSIPN.2018.2869354},
ISSN={2373-776X},
month={June},}

# Requirements

  - MATLAB
  - GSP Toolbox: https://epfl-lts2.github.io/gspbox-html/
  - SPGL1: https://www.cs.ubc.ca/~mpf/spgl1/

GSP Toolbox version 0.7.0 and SPGL1 1.9 provided with the code

# Usage
Launch demo.m to sample and reconstruct a graph signal. Parameters:
```
N = 100;                                % signal length
k = 10;                                 % sparsity
M = 30;                                 % no. of measurements
supp = 1:k;                             % lowpass support
%supp = randperm(N); supp = supp(1:k);  % random sparsity support
graph_name = 'sensor';                  % 'erdos', 'minnesota', 'ring', 'bunny', 'full', 'sensor','scale_free'
supp_known = 1;                         % is the support known?
sigma_noise = 1e-5;                     % standard deviation of noise
```
