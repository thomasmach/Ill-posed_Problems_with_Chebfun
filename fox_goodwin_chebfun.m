function [ke,f,g] = fox_goodwin_chebfun()
% FOX-GOODWIN Test problem: severely ill-posed problem.
%
% [ke,f,g] = fox_goodwin_chebfun()
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
%          1
%  g(s) =  ∫   ke(s,t) f(t) dt
%          0
%
% A first-kind Fredholm integral equation with kernel
%    ke(s,t) = sqrt(s^2 + t^2)
% and with integration intervals
% s \in [0,1],  t \in [0,1].
% The solution is given by
%    f(t) = t.
%
% The right-hand side is
%     g(s) = 1/3 ((1+s^2)^(3/2) - s^3)
%
% This is a model problem which does not satisfy the
% discrete Picard condition for the small singular values.
% The problem was first used by Fox & Goodwin.
%
% Reference: C. T. H. Baker, "The Numerical Treatment of
% Integral Equations", Clarendon Press, Oxford, 1977; p. 665.
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

ke = chebfun2(@(s,t) (s^2+t^2)^(1/2) ,[0 1 0 1],'eps',1e-16,'vectorize');
if (nargout>=2)
  f = chebfun(@(t) t,[0 1],'eps',1e-16,'vectorize');
end
if (nargout>=3)
  g = chebfun(@(s) 1/3*((1+s^2)^(3/2)-s^3) ,[0 1],'eps',1e-16,'vectorize');
end


%%% Local Variables: 
%%% mode:matlab
%%% flyspell-mode:nil
%%% mode:flyspell-prog
%%% ispell-local-dictionary: "american"
%%% End: 
