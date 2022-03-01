function [u12f, si, v12f] = svd4(ke, d)

% INPUT:
%
% ke ......... function handle of ke(s1,s2,t1,t2)
%
% d .......... dimensions of s1, s2, t1, and t2 [s1_lb s1_ub s2_lb s2_ub ... ]
%
% OUTPUT:
%
% ke(s1,s2,t1,t2) = \sum_i u12_i(s1,s2) s_i v12_i(t1,t2)
% 
% 
np = 5000;
tol = 1e-4;



% draw np random points in the domain

s1r = (d(2)-d(1))*rand(np,1) + d(1);
s2r = (d(4)-d(3))*rand(np,1) + d(3);
t1r = (d(6)-d(5))*rand(np,1) + d(5);
t2r = (d(8)-d(7))*rand(np,1) + d(7);


clear function_value_random_points;
for ii=np:-1:1
	function_value_random_points(ii,1) = ke(s1r(ii),s2r(ii),t1r(ii),t2r(ii));
end

ins = ones(np,1);

kf = @(s1,s2,t1,t2) 0;

[ma, ind] = max(abs(function_value_random_points));
mai = ma;
s = [];

kk = 0;
%[function_value_random_points ins]
II = [];

current_function_value_random_points = function_value_random_points;

while (ma>tol*mai)
	kk = kk + 1;
	
	%ind
	
	% compute cross
	u12n = chebfun2(@(s1,s2) ke(s1,s2,t1r(ind),t2r(ind)), [d(1) d(2) d(3) d(4)],'eps',1e-16,'vectorize','splitting','on');
	v12n = chebfun2(@(t1,t2) ke(s1r(ind),s2r(ind),t1,t2), [d(5) d(6) d(7) d(8)],'eps',1e-16,'vectorize','splitting','on');	
	if (kk == 1)
		clear u12;
		clear v12;
		u12{1} = u12n;
		v12{1} = v12n;

	else
		u12{kk} = u12n;
		v12{kk} = v12n;

	end
	
	% cross point
	snew = u12n(s1r(ind),s2r(ind));
	if (kk == 1)
		s = snew;
		u12_discrete = u12n(s1r,s2r); 
 		v12_discrete = v12n(t1r,t2r); 
		
	else
		s = [s, u12n(s1r(II),s2r(II)); transpose(v12n(t1r(II),t2r(II))), snew];
		
		u12_discrete = [u12_discrete, u12n(s1r,s2r)]; 
 		v12_discrete = [v12_discrete, v12n(t1r,t2r)]; 
		
	end
	II = [II ind];
	
	
	% update kf
	kf = @(s1,s2,t1,t2) kff(s1,s2,t1,t2,kk,s,u12,v12);
	

	% update function_value_random_points
	current_function_value_random_points = function_value_random_points - diag(u12_discrete/s*(transpose(v12_discrete)));
	ins(ind) = 0;
	[ma, ind] = max(abs(current_function_value_random_points));
	%[current_function_value_random_points ins]

	check = [current_function_value_random_points(not(ins))];
	if (max(abs(check))>1e-10)
		break
	end
	
end




% qr decomposition of u12 with modified Gram-Schmidt
Ru = zeros(kk,kk);
for ii = 1:kk
	% orthogonalize against previous columns
	for ll = 1:2
		for jj = 1:ii-1
			t = sum2(u12{ii}.*u12{jj});
			if (abs(t)>10*eps)
				u12{ii} = u12{ii} - t*u12{jj};
				Ru(jj,ii) = Ru(jj,ii) + t;
			end
		end
	end
	% normalize
	t = norm(u12{ii});
	if (abs(t-1)>10*eps)
		Ru(ii,ii) = t;
		u12{ii} = u12{ii}/Ru(ii,ii);
	else
		Ru(ii,ii) = 1;
	end
end


% qr decomposition of v with modified Gram-Schmidt
Rv = zeros(kk,kk);
for ii = 1:kk
	% orthogonalize against previous columns
	for ll = 1:2
		for jj = 1:ii-1
			t = sum2(v12{ii}.*v12{jj});
			if (abs(t)>10*eps)
				v12{ii} = v12{ii} - t*v12{jj};
				Rv(jj,ii) = Rv(jj,ii) + t;
			end
		end
	end
	% normalize
	t = norm(v12{ii});
 	if (abs(t-1)>10*eps)
		Rv(ii,ii) = t;
		v12{ii} = v12{ii}/Rv(ii,ii);
	else
		Rv(ii,ii) = 1;
	end
end

%keyboard
[ui,si,vi] = svd(Ru*inv(s)*Rv',0);
%si



%% form u12 * ui
for ii = 1:kk
	t = u12{1}*ui(1,ii);
	for jj=2:kk
		t = t + u12{jj}*ui(jj,ii);
	end
	u12f{ii} = t;
end

for ii = 1:kk
	t = v12{1}*vi(1,ii);
	for jj=2:kk
		t = t + v12{jj}*vi(jj,ii);
	end
	v12f{ii} = t;
end





%%% Local Variables: 
%%% mode:matlab
%%% flyspell-mode:nil
%%% mode:flyspell-prog
%%% ispell-local-dictionary: "american"
%%% End: 
