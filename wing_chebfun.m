function [ke,f,g] = wing_chebfun(t1,t2)
% WING Test problem with a discontinuous solution.
%
% ke           chebfun2
%                the function ke(s,t)
%
% f            chebfun 
%                the function f(t)
% 
% g            chebfun 
%                the function g(s)
%
% discretized outputs
%
% A            matrix of size n x n
%                discretized integral operator 
%
% x            vector of length n
%                discretized exact solution
%
% b            vector of length n
%                discretized right hand side
%
% inputs
% t1, t2       real numbers in (0,1).
%
%          1
%  g(s) =  ∫   ke(s,t) f(t) dt
%          0
% 
%  ke(s,t) = t*exp(-s*t^2)
%         | 1,  t \in [t1,t2],
%  f(t) = {
%         | 0,  else
% 
%  g(s) = (exp(-s*t1^2)-exp(-s*t2^2))/(2*s)
%
% Here, t1 and t2 are constants satisfying t1 < t2.  If they are
% not speficied, the values t1 = 1/3 and t2 = 2/3 are used.
%
% Reference: G. M. Wing, "A Primer on Integral Equations of the
% First Kind", SIAM, 1991; p. 109.
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
%     * Neither the name of the DTU Compute, nor the name of Kent State
%       University, nor the name of the University of Potsdam, nor the names
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
if (nargin<2)
  t1 = 1/3; t2 = 2/3;
else
  if (t1 > t2), error('t1 must be smaller than t2'), end
end

ke = chebfun2(@(s,t) t*exp(-s*t^2),[0 1 0 1],'eps',1e-16,'vectorize');
if (nargout>=2)
  f = chebfun(@(t) ((t1<t)&&(t<t2)),[0 1],'eps',1e-16,'vectorize','splitting','on');
end
if (nargout>=3)
  g = chebfun(@(s) (exp(-s*t1^2)-exp(-s*t2^2))/(2*s),[0 1],'eps',1e-16,'vectorize');
end

