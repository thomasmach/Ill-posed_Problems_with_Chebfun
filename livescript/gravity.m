function [A,b,x,t] = gravity(n,example,a,b,d)
%% GRAVITY Test problem: 1-D gravity surveying model problem
%
% [A,b,x] = gravity(n,example,a,b,d)
%
% Discretization of a 1-D model problem in gravity surveying, in which
% a mass distribution f(t) is located at depth d, while the vertical
% component of the gravity field g(s) is measured at the surface.
%
% The resulting problem is a first-kind Fredholm integral equation
% with kernel
%    K(s,t) = d*(d^2 + (s-t)^2)^(-3/2) .
% The following three examples are implemented (example = 1 is default):
%    1: f(t) = sin(pi*t) + 0.5*sin(2*pi*t),
%    2: f(t) = piecewise linear function,
%    3: f(t) = piecewise constant function.
% The problem is discretized by means of the midpoint quadrature rule
% with n points, leading to the matrix A and the vector x.  Then the
% right-hand side is computed as b = A*x.
%
% The t integration interval is fixed to [0,1], while the s integration
% interval [a,b] can be specified by the user. The default interval is
% [0,1], leading to a symmetric Toeplitz matrix.
%
% The parameter d is the depth at which the magnetic deposit is located,
% and the default value is d = 0.25. The larger the d, the faster the
% decay of the singular values.
%
% Reference: G. M. Wing and J. D. Zahrt, "A Primer on Integral Equations
% of the First Kind", SIAM, Philadelphia, 1991; p. 17.
%
% Per Christian Hansen, IMM, November 18, 2001.
%
% Copyright (c) 2020, Abdulaziz Alqahtani, Lothar Reichel, Thomas Mach
% This file has been modified. It was originally published by
% Copyright (c) 2015, Per Christian Hansen
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%     * Neither the name of the DTU Compute nor the names
%       of its contributors may be used to endorse or promote products derived
%       from this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%

% Initialization.
if (nargin<2), example = 1; end
if (nargin<4), a = 0; b = 1; end
if (nargin<5), d = 0.25; end
if isempty(example), example = 1; end
if isempty(a), a = 0; end
if isempty(b), b = 1; end

% Set up abscissas and matrix.
dt = 1/n;
ds = (b-a)/n;
t = dt*((1:n)' - 0.5);
s = a + ds*((1:n)' - 0.5);
[T,S] = meshgrid(t,s);
A = dt*d*ones(n,n)./(d^2 + (S-T).^2).^(3/2);

% Set up solution vector and right-hand side.
nt = round(n/3);
nn = round(n*7/8);
x = ones(n,1);
switch example
case 1
   x = sin(pi*t) + 0.5*sin(2*pi*t);
case 2
   x(1:nt)    = (2/nt)*(1:nt)';
   x(nt+1:nn) = ((2*nn-nt) - (nt+1:nn)')/(nn-nt);
   x(nn+1:n)  = (n - (nn+1:n)')/(n-nn);
case 3
   x(1:nt) = 2*ones(nt,1);
otherwise
   error('Illegal value of example')
end
b = A*x;