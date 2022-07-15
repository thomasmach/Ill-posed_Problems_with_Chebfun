function [A,b,x,t] = wing(n,t1,t2)
%% WING Test problem with a discontinuous solution.
%
% [A,b,x] = wing(n,t1,t2)
%
% Discretization of a first kind Fredholm integral eqaution with
% kernel K and right-hand side g given by
%    K(s,t) = t*exp(-s*t^2)                       0 < s,t < 1
%    g(s)   = (exp(-s*t1^2) - exp(-s*t2^2)/(2*s)  0 < s   < 1
% and with the solution f given by
%    f(t) = | 1  for  t1 < t < t2
%           | 0  elsewhere.
%
% Here, t1 and t2 are constants satisfying t1 < t2.  If they are
% not speficied, the values t1 = 1/3 and t2 = 2/3 are used.
%
% Reference: G. M. Wing, "A Primer on Integral Equations of the
% First Kind", SIAM, 1991; p. 109.
%
% Discretized by Galerkin method with orthonormal box functions;
% both integrations are done by the midpoint rule.
%
% Per Christian Hansen, IMM, 09/17/92.
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
if (nargin==1)
  t1 = 1/3; t2 = 2/3;
else
  if (t1 > t2), error('t1 must be smaller than t2'), end
end
A = zeros(n,n); h = 1/n;

% Set up matrix.
sti = ((1:n)-0.5)*h;
for i=1:n
  A(i,:) = h*sti.*exp(-sti(i)*sti.^2);
end

% Set up right-hand side.
if (nargout > 1)
  b = sqrt(h)*0.5*(exp(-sti*t1^2)' - exp(-sti*t2^2)')./sti';
end

% Set up solution.
if (nargout>=3)
  I = find(t1 < sti & sti < t2);
  x = zeros(n,1); x(I) = sqrt(h)*ones(length(I),1);
end

% provide the used discretization.
if (nargout==4)
	t = (.5:n-.5)'*h;
end