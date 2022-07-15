function [err] = two_errlambda1(lambda,sigma,gnoise,psi,rk,eta,e,dd)
%% Error function for fminbnd optimization -- continuous rhs
% The purpose of this function is to evaluate
% || Ax - g || + lambda || x || - delta
%
% [err] = twoerrlambda1(lambda,sigma,gnoise,phi,rk,eta,e,dd)
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

beta = dd.*sigma./(sigma.^2+lambda^2);  

kex = beta(1)*sigma(1)*psi{1};
for ll = 2:rk
	kex = kex + beta(ll)*sigma(ll)*psi{ll};
end
err = abs(norm(kex - gnoise)^2-eta^2*e^2);



%%% Local Variables: 
%%% mode:matlab
%%% flyspell-mode:nil
%%% mode:flyspell-prog
%%% ispell-local-dictionary: "american"
%%% End: 
