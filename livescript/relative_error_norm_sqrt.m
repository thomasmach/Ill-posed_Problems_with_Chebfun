function r = relative_error_norm_sqrt(f,t,x,opt)
%% Function to compute the relative error norm 
% of a function defined by the discrete vector x at timepoint t and the chebfun
% f.
%
% If opt is true, then linear interpolation is used. Otherwise constant
% functions are assumed.

debug = false;


if (nargin<4)
	opt = false;
end

n = length(t);
m = length(x);
assert(n==m);
ftt = f.domain;
ft = [ftt(1) ftt(end)];

assert(prod(t>=ft(1))==1);
assert(prod(t<=ft(2))==1);

if (debug)
	figure
	hold on
	plot(f)
end

ii = 1;
if (opt)
	
	g = chebfun(@(y) (x(1)+(y-t(1))/(t(2)-t(1))*(x(2)-x(1)))/sqrt(0.5*(t(1)+t(2))-ft(1))-f(y), [ft(1), t(2)],'splitting','on');
	if (debug)
		h = chebfun(@(y) (x(1)+(y-t(1))/(t(2)-t(1))*(x(2)-x(1)))/sqrt(0.5*(t(1)+t(2))-ft(1)), [ft(1), t(2)],'splitting','on');
		plot(g)
		plot(h)
	end
	
else
	g = chebfun(@(y) x(1)/sqrt(0.5*(t(1)+t(2))-ft(1))-f(y), [ft(1), 0.5*(t(1)+t(2))],'splitting','on');
	if (debug)
		h = chebfun(@(y) x(1)/sqrt(0.5*(t(1)+t(2))-ft(1)), [ft(1), 0.5*(t(1)+t(2))],'splitting','on');
		plot(g)
		plot(h)
	end
end
r = norm(g)^2;

for ii=2:n-2
	
	if (opt)
		g = chebfun(@(y) (x(ii) + (y-t(ii))/(t(ii+1)-t(ii))*(x(ii+1)-x(ii)))/sqrt(0.5*(t(ii+1)-t(ii-1)))-f(y), [t(ii), t(ii+1)],'splitting','on');
		if (debug)
			h = chebfun(@(y) (x(ii) + (y-t(ii))/(t(ii+1)-t(ii))*(x(ii+1)-x(ii)))/sqrt(0.5*(t(ii+1)-t(ii-1))), [t(ii),t(ii+1)],'splitting','on');
			plot(g)
			plot(h)
		end
	else
		g = chebfun(@(y) x(ii)/sqrt(0.5*(t(ii+1)-t(ii-1)))-f(y), [0.5*(t(ii-1)+t(ii)), 0.5*(t(ii)+t(ii+1))],'splitting','on');
		if (debug)
			h = chebfun(@(y) x(ii)/sqrt(0.5*(t(ii+1)-t(ii-1))), [0.5*(t(ii-1)+t(ii)), 0.5*(t(ii)+t(ii+1))],'splitting','on');
			plot(g)
			plot(h)
		end
	end
	r = r+norm(g)^2;
	
end

ii = n;
if (opt)

	g = chebfun(@(y) (x(n-1) + (y-t(n-1))/(t(n)-t(n-1))*(x(n)-x(n-1)))/sqrt(ft(2)-0.5*(t(n-1)+t(n)))-f(y), [t(n-1), ft(2)],'splitting','on');
	if (debug)
		h = chebfun(@(y) (x(n-1) + (y-t(n-1))/(t(n)-t(n-1))*(x(n)-x(n-1)))/sqrt(ft(2)-0.5*(t(n-1)+t(n))), [t(n-1), ft(2)],'splitting','on');
		plot(g)
		plot(h)
	end	
	
else

	ii = n - 1;	
	if (n>2)
		g = chebfun(@(y) x(ii)/sqrt(0.5*(t(ii+1)-t(ii-1)))-f(y), [0.5*(t(ii-1)+t(ii)), 0.5*(t(ii)+t(ii+1))],'splitting','on');
		if (debug)
			h = chebfun(@(y) x(ii)/sqrt(0.5*(t(ii+1)-t(ii-1))), [0.5*(t(ii-1)+t(ii)), 0.5*(t(ii)+t(ii+1))],'splitting','on');
			plot(g)
			plot(h)
		end
		r = r+norm(g)^2;
	end
	
	ii = n;
	g = chebfun(@(y) x(n)/sqrt(ft(2)-0.5*(t(n-1)+t(n)))-f(y), [0.5*(t(n-1)+t(n)), ft(2)],'splitting','on');
	if (debug)
		h = chebfun(@(y) x(n)/sqrt(ft(2)-0.5*(t(n-1)+t(n))), [0.5*(t(n-1)+t(n)), ft(2)],'splitting','on');
		plot(g)
		plot(h)
	end

end

r = r+norm(g)^2;
r = sqrt(r) / norm(f);

if (debug)
	keyboard
	hold off
end



%%% Local Variables: 
%%% flyspell-mode:nil
%%% mode:flyspell-prog
%%% mode:matlab
%%% ispell-local-dictionary: "american"
%%% End: 
