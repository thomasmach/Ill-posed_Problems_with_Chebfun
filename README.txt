********************************************************************************
*                                                                              * 
* MIT License                                                                  *
*                                                                              * 
* Copyright (c) 2022  Abdulaziz Alqahtani, Thomas Mach, Lothar Reichel         *
*                                                                              * 
* Permission is hereby granted, free of charge, to any person obtaining a copy *
* of this software and associated documentation files (the "Software"), to deal*
* in the Software without restriction, including without limitation the rights *
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell    *
* copies of the Software, and to permit persons to whom the Software is        *
* furnished to do so, subject to the following conditions:                     *
*                                                                              * 
* The above copyright notice and this permission notice shall be included in   *
* all copies or substantial portions of the Software.                          *
*                                                                              * 
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR   *
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,     *
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE  *
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER       *
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,*
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE*
* SOFTWARE.                                                                    *
*                                                                              * 
********************************************************************************
*                                                                              * 
* Some files are modified files originally published by Per Christian Hansen   *
* (DTU) under a 3-Clause BSD License. The license for these files is stated at *
* the begin of the relevant files.                                             *
*                                                                              * 
********************************************************************************

 AUTHORS:

     ABDULAZZIZ ALQAHTANI
		 KINK KHALID UNIVERSITY - SAUDI ARABIA
		 KENT STATE UNIVERSITY - USA
		 E-MAIL: aalqah11@gmail.com

		 THOMAS MACH
		 UNIVERSITY OF POTSDAM - GERMANY
		 E-MAIL: mach@uni-potsdam.de

		 LOTHAR REICHEL
		 KENT STATE UNIVERSITY - USA
		 E-MAIL: reichel@math.kent.edu

 REFERENCE:

  -  Solution of Ill-Posed Problems with Chebfun
		 NUMERICAL ALGORITHMS, TBA

 SOFTWARE REVISION:

     V1.0, May 2022

 SOFTWARE LANGUAGE:

     Matlab

 REQUIREMENTS:

     Chebfun 5.7.0 or newer (available from https://github.com/chebfun/chebfun)

********************************************************************************

This folder contains the following files:

* README.txt

* README.md
  landing page for github users

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
