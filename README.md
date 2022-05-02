# Solution of Ill-posed Problems with Chebfun

This repository contains MATLAB livescripts and files corresponding to the paper
“Solution of Ill-posed Problems with Chebfun” by Abdulaziz Alqahtani, Thomas
Mach, Lothar Reichel, Numerical Algorithms, TBA.

The code was tested with MATLAB R2020b and the latest release of Chebfun
(v.5.7.0). The Chebfun package is [available on github under https://github.com/chebfun/chebfun](https://github.com/chebfun/chebfun).

The examples are mostly modified version of [Regularization Tools Version 4.1](http://www.imm.dtu.dk/~pcha/Regutools/).
One example is inspired by 

Gazzola, S., Hansen, P. C., & Nagy, J. G. (2019). IR Tools: a MATLAB package
of iterative regularization methods and large-scale test problems. Numerical
Algorithms, 81(3), 773-811.


This folder contains the following files:

* LICENSE
	license file

* baart.m, fox_goodwin.m, gravity.m, shaw.m, wing.m
	definition of examples from Regularization Tools by Per Christian Hansen
	the files have been modified to provide the used discretization

* baart_chebfun.m, fox_goodwin_chebfun.m gravity_chebfun.m, shaw_chebfun.m
  wing_chebfun.m
	continuous version of the examples above using Chebfun

* Solution_of_Ill_posed_Problems_with_Chebfun_1d/2d_live.mlx
	Matlab livescripts demonstrating how the examples can be solved

* svd4.m
	singular value decomposition for four dimensional functions using
	adaptive cross approximation
	
* errlambda1/2b/3.m, twod_errlambda1.m
	support function
	different version of objective functions for the different examples

* relative_error_norm.m/relative_error_norm_sqrt.m
	support function

* weighted_euclidean_norm.m
	support function
