function [ke,f,g] = baart_chebfun(f0)
%% BAART Test problem: Fredholm integral equation of the first kind.
%
% [ke,f,g] = baart()
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
% f0           complex
%                constant added to original f(t), default 0
%
%         \pi
%  g(s) =  ∫   ke(s,t) f(t) dt
%          0
%
% A first-kind Fredholm integral equation with kernel
%    ke(s,t) = exp(s*cos(t))
% and with integration intervals
%  s \in [0,pi/2],  t \in [0,pi].
% The solution is given by
%    f(t) = f0 + sin(t).
%
% If f0==0 then the right hand side is 
%    g(s) = 2*sinh(s)/s.
%
% Reference: M. L. Baart, "The use of auto-correlation for pseudo-
% rank determination in noisy ill-conditioned linear least-squares
% problems", IMA J. Numer. Anal. 2 (1982), 241-247. Example 4.2
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
  f0 = 0;
end

ke = chebfun2(@(s,t) exp(s*cos(t)) ,[0 pi/2 0 pi],'eps',1e-16,'vectorize');
if (nargout>=2)
  f = chebfun(@(t) f0+sin(t),[0 pi],'eps',1e-16,'vectorize');
end
if (nargout>=3)
  if (f0==0)
    g = chebfun(@(s) 2*sinh(s)/s,[0 pi/2],'eps',1e-16,'vectorize');
  else
    g = sum(transpose(ke)*f,2);
  end
end

