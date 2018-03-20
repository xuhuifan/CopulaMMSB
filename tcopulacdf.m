function c_val = tcopulacdf( u,v,ctheta, ctheta_2 )
% Gaussian copula cdf
%   u, v :  two values
%   ctheta: copula function parameter
ro = [1 ctheta;ctheta 1];
u1 = reshape(u, [], 1);
v1 = reshape(v, [], 1);
c_val = mvtcdf(tinv([u1 v1],ctheta_2),ro,ctheta_2);
c_val = reshape(c_val, size(u));
end

