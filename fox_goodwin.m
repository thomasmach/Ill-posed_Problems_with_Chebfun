function [A,b,x,t] = fox_goodwin(n)
%% FOX-GOODWIN (aka FOXGOOD) Test problem: severely ill-posed problem.
%
% [A,b,x] = fox_goodwin(n)
%
% This is a model problem which does not satisfy the
% discrete Picard condition for the small singular values.
% The problem was first used by Fox & Goodwin.
%
% Reference: C. T. H. Baker, "The Numerical Treatment of
% Integral Equations", Clarendon Press, Oxford, 1977; p. 665.
%
% Discretized by simple quadrature (midpoint rule).
%
% Per Christian Hansen, IMM, 03/16/93.
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
h = 1/n; t = h*((1:n)' - 0.5);

A = h*sqrt((t.^2)*ones(1,n) + ones(n,1)*(t.^2)');
x = t; b = ((1+t.^2).^1.5 - t.^3)/3;

