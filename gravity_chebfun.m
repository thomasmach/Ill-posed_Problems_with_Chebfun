function [ke,f,g] = gravity_chebfun(example,d)
% GRAVITY Test problem: 1-D gravity surveying model problem
% [ke,f,g] = gravity(n,example,a,b,d)
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
% example      integer \in {1,2,3}
%                example no., default 1
%
% d            complex
%                parameter, default 0.25 
%                (larger d results in a smaller rank of ke)
%
%          1
%  g(s) =  ∫   ke(s,t) f(t) dt
%          0
%
% 1-D model problem in gravity surveying, in which a mass distribution f(t) is
% located at depth d, while the vertical component of the gravity field g(s) is
% measured at the surface.  
%
% The resulting problem is a first-kind Fredholm integral equation with kernel
%    ke(s,t) = d*(d^2 + (s-t)^2)^(-3/2) .
% The following three examples are implemented (example = 1 is default):
%    example = 1: f(t) = sin(pi*t) + 0.5*sin(2*pi*t),
%    example = 2: f(t) = piecewise linear function,
%    example = 3: f(t) = piecewise constant function.
%
% Both integration interval are [0,1].
%
% The parameter d is the depth at which the magnetic deposit is located,
% and the default value is d = 0.25. The larger the d, the faster the
% decay of the singular values.
%
% Reference: G. M. Wing (with the assistance of J. D. Zahrt), "A Primer on
% Integral Equations of the First Kind", SIAM, Philadelphia, 1991; p. 17.
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

if (nargin<1) 
  example = 1; 
end
if (nargin<2) 
  d = 0.25; 
end
if (isempty(example))
  example = 1; 
end

ke = chebfun2(@(s,t) d*(d^2 + (s-t)^2)^(-3/2),[0 1 0 1],'eps',1e-16,'vectorize','splitting','on');

if (nargout>=2)
  switch (example)
    case 1
      f = chebfun(@(t) sin(pi*t) + 0.5*sin(2*pi*t),[0 1],'eps',1e-16,'vectorize','splitting','on');
      
    case 2
      f = chebfun(@(t) 6*t*(t<=1/3) + (2-(t-1/3)*24/13)*((t<7/8)&&(t>1/3)) + 8*(1-t)*(t>=7/8) ,[0 1],'eps',1e-16,'vectorize','splitting','on');
      
    case 3
      f = chebfun(@(t) 1+(t<1/3),[0 1],'eps',1e-16,'vectorize','splitting','on');
    otherwise
      error('Illegal value of example')
  end
end
if (nargout>=3)
  g = sum(ke*f,2);
end

%%% Local Variables: 
%%% mode:matlab
%%% flyspell-mode:nil
%%% mode:flyspell-prog
%%% ispell-local-dictionary: "american"
%%% End: 
