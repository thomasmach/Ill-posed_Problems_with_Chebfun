%%% Error function for fminbnd optimization -- continuous rhs
% The purpose of this function is to evaluate
% || Ax - g || + lambda || x || - delta
%
% XXXX add discription what parameterchoice exactly does XXXX
% 
function err=errlambda3(parameterchoice,lambda,rhs,phi,gnoise,psi,ke,rk,eta,e,dd,E,D)

DD = E + lambda^2*D;
beta3=DD\rhs;

x = psi(:,1:rk)*beta3(1:rk,1); 

switch (parameterchoice)

	case {1}
		err = abs(norm(sum(ke'*x,2)-gnoise)^2-eta^2*e^2);
		
	case {2}
		err = (lambda^(4)*dd.^2)./(sigma.^2+lambda^2).^2;
		err = abs(sum(err)-(eta*e)^2);

	case {3}
		err = lambda^6*dd'*((sigma.^2+lambda^2).^(-3).*dd);
		err=abs((err)-(eta*e)^2);
end


end


%%% Local Variables: 
%%% mode:matlab
%%% flyspell-mode:nil
%%% mode:flyspell-prog
%%% ispell-local-dictionary: "american"
%%% End: 
