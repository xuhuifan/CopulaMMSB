function c_val = gcopulacdf( u,v,ctheta )
% Gaussian copula cdf
%   u, v :  two values
%   ctheta: copula function parameter
ro = [1 ctheta;ctheta 1];
u1 = reshape(u, [], 1);
v1 = reshape(v, [], 1);
c_val = mvncdf(norminv([u1 v1]),zeros(1,2),ro);
c_val = reshape(c_val, size(u));
end



