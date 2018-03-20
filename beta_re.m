function betas = beta_re(cpi, cu_betas, numClass, dataNum, alphas)

new_betas = dirrnd(2*ones(1, numClass), 1);
c_value = sum((alphas*cu_betas).*sum(log(cpi)))-dataNum*sum(gammaln(alphas*cu_betas));

q_value = sum((alphas*new_betas).*sum(log(cpi)))-dataNum*sum(gammaln(alphas*new_betas));

if rand < (q_value/c_value)
    betas = new_betas;
else
    betas = cu_betas;
end





