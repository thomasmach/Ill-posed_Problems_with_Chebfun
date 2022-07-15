function [err] = errlambda1(lambda,sigma,gnoise,psi,ke,rk,eta,e,dd)
%% Error function for fminbnd optimization -- continuous rhs
% The purpose of this function is to evaluate
% | || Ax - g || + lambda || x || - eta * e |
% This is used as a cost function for the fminbnd optimization to find a suitable
% lambda.
%
% [err] = errlambda1(lambda,sigma,gnoise,psi,ke,rk,eta,e,dd)
%
% INPUT:
%
% lambda ............. double
%                        regularization parameter to be optimized
%
% sigma .............. vector
%                        singular values
%
% gnoise ............. chebfun
%                        noisy right-hand side
%
% psi ................ vector of chebfun
%                        left singular functions
% 
% ke ................. chebfun2
%                        kernel function
%
% rk ................. integer
%                        rank, include the first rk singular values/vectors
%
% eta ................ double
%                        additional factor in the discrepancy principle, 
%                        choose between 1 and 5
%
% e .................. double
%                        estimation of the norm of the noise contained in gnoise
%
% dd ................. vector
%                        inner product of gnoise with the right singular functions phi
%
% OUTPUT:
%
% err ................ double 
%                        value of the cost function
% 

% computng the solution
beta = dd.*sigma./(sigma.^2+lambda^2);  
x = psi(:,1:rk)*beta(1:rk,1); 

% evaluating the cost function
err = abs(norm(sum(transpose(ke)*x,2)-gnoise)^2-eta^2*e^2);



%%% Local Variables: 
%%% mode:matlab
%%% flyspell-mode:nil
%%% mode:flyspell-prog
%%% ispell-local-dictionary: "american"
%%% End: 
