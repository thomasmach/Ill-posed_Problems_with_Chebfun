function [ke,f,g] = shaw_chebfun()
%% SHAW Test problem: one-dimensional image restoration model.
%
% [ke,f,g] = shaw_chebfun()
%
% INPUT: NONE
%
% OUTPUT:
%
% ke ......... chebfun2
%                the function ke(s,t)
%
% f .......... chebfun 
%                the function f(t)
% 
% g .......... chebfun 
%                the function g(s)
%
%
%         \pi/2
%  g(s) =  ∫   ke(s,t) f(t) dt
%        -\pi/2
%
% A first-kind Fredholm integral equation with kernel
%    ke(s,t) = (cos(s)+cos(t))*(sin(u)/u)^2, with 
%          u = pi*(sin(s)+sin(t))
% and with integration intervals
% s \in [-pi/2,pi/2],  t \in [-pi/2,pi/2].
% The solution is given by
%    f(t) = a1*exp(-c1*(t - t1)^2) + a2*exp(-c2*(t - t2)^2), with
%      a1 = 2; c1 = 6; t1 =  .8;
%      a2 = 1; c2 = 2; t2 = -.5;
%
% The right-hand side is computed based on ke and f.
%
% Reference: C. B. Shaw, Jr., "Improvements of the resolution of
% an instrument by numerical solution of an integral equation",
% J. Math. Anal. Appl. 37 (1972), 83-112. In particular pp. 97 and 98.
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


u = chebfun2(@(s,t) pi()*(sin(s)+sin(t)),[-pi()/2 pi()/2 -pi()/2 ...
      pi()/2],'eps',1e-16,'splitting','on');

q = chebfun(@(x) (sin(x)/x).^2, [-2,2],'vectorize');
qu = chebfun2(@(s,t) q(u(s,t)), ...
    [-pi()/2 pi()/2 -pi()/2 pi()/2],'eps',1e-16,'splitting','on');

cs = chebfun2(@(s,t) (cos(s)+cos(t)),[-pi()/2 pi()/2 -pi()/2 pi()/2],'eps',1e-16);

ke = cs.*qu;

a1 = 2;
a2 = 1;
c1 = 6;
c2 = 2;
t1 =  .8;
t2 = -.5;

if (nargout>=2)
  f = chebfun(@(t) a1*exp(-c1*(t-t1)^2)+a2*exp(-c2*(t-t2)^2),[-pi()/2 pi()/2],'eps',1e-16,'vectorize');
end
if (nargout>=3)
  g = sum(transpose(ke)*f,2);
end

%%% Local Variables: 
%%% mode:matlab
%%% flyspell-mode:nil
%%% mode:flyspell-prog
%%% ispell-local-dictionary: "american"
%%% End: 

