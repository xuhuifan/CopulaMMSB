function [cthetas cthetas_2]= hyper_para( dataNum, cpi, seLabel, reLabel, ctheta, ctheta_2, de_num)
% sampling the hyper-parameters \gamma, \alpha

% input:
%  m_val: K-length vector denoting each component's table number
%  nums: N-length vector denoting each restaurant's dish number
% dataNum: the data points number
% cpi : the \pi's value
% seLabel : \{s_{ij}\}_{n\times n}
% reLabel : \{r_{ij}\}_{n\times n}
% ctheta : the currrent copula function parameter

% output:
% alphas: the required alpha value
% gammas: the required gamma value
% cthetas: the new copula function parameter theta's value
%
%
new_theta = gamrnd(1, 1)+1;
theta_ratio = theta_ar_1( cpi, seLabel, reLabel, ctheta, new_theta, dataNum, de_num);
%  fprintf('theta_ratio is %f\n', theta_ratio);
if rand < theta_ratio
    cthetas = new_theta;
else
    cthetas = ctheta;
end

new_theta = gamrnd(1, 1)+1;
theta_ratio = theta_ar_2( cpi, seLabel, reLabel, ctheta_2, new_theta, dataNum, de_num);
%  fprintf('theta_ratio is %f\n', theta_ratio);
if rand < theta_ratio
    cthetas_2 = new_theta;
else
    cthetas_2 = ctheta_2;
end

end

