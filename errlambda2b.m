function [err2] = errlambda2b(A,bnoise,lambda,s,foriercoeff,V,eta,e)
%% Error function for fminbnd optimization -- discrete rhs
% The purpose of this function is to evaluate
% || Ax - g || + lambda || x || - delta
%

x_lambda = V*(s.*foriercoeff./(s.^2+lambda^2));
err2 = abs(weighted_euclidean_norm(A*x_lambda-bnoise)^2-eta^2*e^2);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% beta2 = sigma./(sigma.^2+lambda^2);
% beta = zeros(rk,1);
% for ii=1:rk
%     beta(ii) = sum(phi(sd,ii).*gnoise,1)*beta2(ii)*h;
% end
% x = psi(:,1:rk)*beta(1:rk,1);
% ga = sum(transpose(ke)*x,2);
% %err = abs(norm(ga(sd)-gnoise)-eta*delta)/delta; 
% err = abs(norm(ga(sd)-gnoise)-eta*e);



%%% Local Variables: 
%%% mode:matlab
%%% flyspell-mode:nil
%%% mode:flyspell-prog
%%% ispell-local-dictionary: "american"
%%% End: 
