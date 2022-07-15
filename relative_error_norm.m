function [r] = relative_error_norm(f,t,x,opt)
%% Function to compute the relative error norm 
% of a function defined by the discrete vector x at timepoints t and the chebfun
% f.
%
% [r] = relative_error_norm(f,t,x,opt)
%
% INPUT:
%
% f .................. chebfun
%                        function f
%
% t .................. vector
%                        discretization points for x
%
% x .................. vector
%                        vector x with which f is compared
%
% opt ................ boolean
%                        true ............. linear interpolation
%                        false (default) .. piecewise constant
%
% OUTPUT:
%
% r .................. double 
%                        relative normwise error
%

% default option
if (nargin<4)
	opt = false;
end

% setup
n = length(t);
m = length(x);
assert(n==m);
ftt = f.domain;
ft = [ftt(1) ftt(end)];

% check that the timepoints t are inside the domain of the chebfun f
assert(prod(t>=ft(1))==1);
assert(prod(t<=ft(2))==1);

% first interval
ii = 1;
if (opt)
	g = chebfun(@(y) (x(1)+(y-t(1))/(t(2)-t(1))*(x(2)-x(1)))-f(y), [ft(1), t(2)],'splitting','on');
else
	g = chebfun(@(y) x(1)-f(y), [ft(1), 0.5*(t(1)+t(2))],'splitting','on');
end
r = norm(g)^2;

% middle intervals
for ii=2:n-2
	if (opt)
		g = chebfun(@(y) (x(ii) + (y-t(ii))/(t(ii+1)-t(ii))*(x(ii+1)-x(ii)))-f(y), [t(ii), t(ii+1)],'splitting','on');
	else
		g = chebfun(@(y) x(ii)-f(y), [0.5*(t(ii-1)+t(ii)), 0.5*(t(ii)+t(ii+1))],'splitting','on');
	end
	r = r+norm(g)^2;
end

% final interval
ii = n;
if (opt)
	g = chebfun(@(y) (x(n-1) + (y-t(n-1))/(t(n)-t(n-1))*(x(n)-x(n-1)))-f(y), [t(n-1), ft(2)],'splitting','on');	
else
	ii = n - 1;	
	if (n>2)
		g = chebfun(@(y) x(ii)-f(y), [0.5*(t(ii-1)+t(ii)), 0.5*(t(ii)+t(ii+1))],'splitting','on');
		r = r+norm(g)^2;
	end
	
	ii = n;
	g = chebfun(@(y) x(n)-f(y), [0.5*(t(n-1)+t(n)), ft(2)],'splitting','on');

end
r = r+norm(g)^2;
r = sqrt(r) / norm(f);



%%% Local Variables: 
%%% flyspell-mode:nil
%%% mode:flyspell-prog
%%% mode:matlab
%%% ispell-local-dictionary: "american"
%%% End: 
