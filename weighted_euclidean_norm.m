function [wn] = weighted_euclidean_norm(x)
%% Weighted Euclidean norm of a vector x\in\Rn
% ||x||/sqrt(n)
%
% [wn] = weighted_euclidean_norm(x)
%
% INPUT:
%
% x .................. vector
%                        vector x of length n
%
% OUTPUT:
%
% wn ................. double 
%                        value ||x||_2/sqrt(n)
% 

n = length(x);
wn = sqrt(1/n) * norm(x,2);
