function alphas = randconparam(alpha,numdata,numclass,aa,bb,numiter);

% Modification of Escobar and West.  Works for multiple groups of data.
% numdata, numclass are row vectors, one element per group.

if nargin == 5
  numiter = 1;
end

totalclass = sum(numclass);
num = length(numdata);
alphas = zeros((1+numiter), 1);
alphas(1) = alpha;
for ii = 1:numiter
  xx = randbeta((alphas(ii)+1)*ones(1,num),numdata);

  zz = rand(1,num).*(alphas(ii)+numdata)<numdata;

  gammaa = aa + totalclass - sum(zz);
  gammab = bb - sum(log(xx));
  alphas(ii+1) = randgamma(gammaa)./gammab;

end
alphas = alphas(2:end);

