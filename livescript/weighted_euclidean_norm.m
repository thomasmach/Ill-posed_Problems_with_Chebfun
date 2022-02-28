function wn = weighted_euclidean_norm(x)
n = length(x);
wn = sqrt(1/n) * norm(x,2);
end