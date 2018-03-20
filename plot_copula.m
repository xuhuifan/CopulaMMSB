N = 10000; 
u = copularnd('Gumbel', 3.5, N);
u = u';
plot(u(1,:),u(2,:), '.');


