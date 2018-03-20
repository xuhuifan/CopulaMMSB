function c_val = gumbelcdf( u,v,ctheta )
% Frank copula cdf
%   u, v :  two values
%   ctheta: copula function parameter
c_val = exp(-((-log(u)).^ctheta+(-log(v)).^ctheta).^(1/ctheta));
end