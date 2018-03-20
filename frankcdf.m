function c_val = frankcdf( u,v,ctheta )
% Frank copula cdf
%   u, v :  two values
%   ctheta: copula function parameter
c_val = -log(1+(exp(-ctheta*u)-1).*(exp(-ctheta*v)-1)/(exp(-ctheta)-1))/ctheta;
end