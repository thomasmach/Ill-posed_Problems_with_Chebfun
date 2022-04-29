%%% Error function for fminbnd optimization -- continuous rhs
% The purpose of this function is to evaluate
% || Ax - g || + lambda || x || - delta
%
% XXXX add discription what parameterchoice exactly does XXXX
% 
function err=twoerrlambda1(parameterchoice,lambda,sigma,gnoise,phi,rk,eta,e,dd)

%keyboard
beta = dd.*sigma./(sigma.^2+lambda^2);  
%x = psi{1}*beta(1);
%for ll = 2:rk
%	x = x + psi{ll}*beta(ll); 
%end

switch (parameterchoice)

	case {1}
		kex = beta(1)*sigma(1)*phi{1};
		for ll = 2:rk
			kex = kex + beta(ll)*sigma(ll)*phi{ll};
		end
		%err = abs(sum2((kex - gnoise).^2)-eta^2*e^2);
         err = abs(norm(kex - gnoise)^2-eta^2*e^2);
% 	case {2}
% 		err = (lambda^(4)*dd.^2)./(sigma.^2+lambda^2).^2;
% 		err = abs(sum(err)-(eta*e)^2);
% 
% 	case {3}
% 		% err=lambda^(4)*dd.^2.*sigma.^2./(sigma.^2+lambda^2).^2;
% 		% err=abs(sum(err)-(eta*e)^2);
% 		err = lambda^6*dd'*((sigma.^2+lambda^2).^(-3).*dd);
% 		
% 		err=abs((err)-(eta*e)^2);
end


end


%%% Local Variables: 
%%% mode:matlab
%%% flyspell-mode:nil
%%% mode:flyspell-prog
%%% ispell-local-dictionary: "american"
%%% End: 
