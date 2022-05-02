function err=errlambda1(parameterchoice,lambda,sigma,gnoise,psi,ke,rk,eta,e,dd)
%% Error function for fminbnd optimization -- continuous rhs
% The purpose of this function is to evaluate
% || Ax - g || + lambda || x || - delta
%
% XXXX add discription what parameterchoice exactly does XXXX
% 
 
beta = dd.*sigma./(sigma.^2+lambda^2);  
x = psi(:,1:rk)*beta(1:rk,1); 

switch (parameterchoice)

	case {1}
		err = abs(norm(sum(transpose(ke)*x,2)-gnoise)^2-eta^2*e^2);
		
	case {2}
		err = (lambda^(4)*dd.^2)./(sigma.^2+lambda^2).^2;
		err = abs(sum(err)-(eta*e)^2);

	case {3}
		err = lambda^6*dd'*((sigma.^2+lambda^2).^(-3).*dd);
		err=abs((err)-(eta*e)^2);

end



%%% Local Variables: 
%%% mode:matlab
%%% flyspell-mode:nil
%%% mode:flyspell-prog
%%% ispell-local-dictionary: "american"
%%% End: 
