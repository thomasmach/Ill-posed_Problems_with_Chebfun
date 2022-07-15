function [err2] = errlambda2b(A,bnoise,lambda,s,foriercoeff,V,eta,e)
%% Error function for fminbnd optimization -- discrete rhs
% The purpose of this function is to evaluate
% The purpose of this function is to evaluate
% | || Ax - g || + lambda || x || - eta * e |
% This is used as a cost function for the fminbnd optimization to find a suitable
% lambda.
%
% [err2] = errlambda2b(A,bnoise,lambda,s,foriercoeff,V,eta,e)
%
% INPUT:
%
% A .................. matrix
%                        coefficient matrix
% 
% bnoise ............. vector
%                        noisy right-hand side
%
% lambda ............. double
%                        regularization parameter to be optimized
%
% s .................. vector
%                        singular values
%
% fouriercoeff ....... vector
%                        inner product of bnoise with the right singular vectors
%
% V .................. vector
%                        left singular functions
% 
% eta ................ double
%                        additional factor in the discrepancy principle, 
%                        choose between 1 and 5
%
% e .................. double
%                        estimation of the norm of the noise contained in gnoise
%
% OUTPUT:
%
% err2 ............... double 
%                        value of the cost function
 

% computing the solution
x_lambda = V*(s.*foriercoeff./(s.^2+lambda^2));

% evaluating the cost function
err2 = abs(weighted_euclidean_norm(A*x_lambda-bnoise)^2-eta^2*e^2);




%%% Local Variables: 
%%% mode:matlab
%%% flyspell-mode:nil
%%% mode:flyspell-prog
%%% ispell-local-dictionary: "american"
%%% End: 
