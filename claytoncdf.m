function c_val = claytoncdf( u,v,ctheta )
% Clayton copula cdf 
%   u, v :  two values
%   ctheta: copula function parameter
c_val = (u.^(-ctheta)+v.^(-ctheta)-1).^(-1/ctheta);
end

