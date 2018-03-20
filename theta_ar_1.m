function theta_ratio = theta_ar( cpi, seLabel, reLabel, ctheta, new_theta, dataNum, de_num)
% calcluate the acceptance ratio
% cpi : the current whole distribution
% seLabel: \{s_{ij}\}_{i,j}
% reLabel: \{r_{ij}\}_{i,j}
% ctheta: copula function parameter
% new_theta: the proposed new copula function parameter
% dataNum : data size

cu_pi = [zeros(dataNum, 1) cumsum(cpi, 2)]; % the cumulated value \sum_{d=1}^k \pi_{id}
de_len = length(de_num);
recu22 = zeros(de_len, de_len);
recv22 = zeros(de_len, de_len);
recu11 = zeros(de_len, de_len);
recv11 = zeros(de_len, de_len);
% for i = 1:dataNum
%     cu22(i,:) = cu_pi(i,(seLabel(i,:)+1))';
%     cv22(i,:) = diag(cu_pi(1:dataNum,(reLabel(i,:)+1)));
%     cu11(i,:) = cu_pi(i,(seLabel(i,:)))';
%     cv11(i,:) = diag(cu_pi(1:dataNum,(reLabel(i,:))));
% end
% qc = (indecdf(cu22, cv22)-indecdf(cu11, cv22)+indecdf(cu11, cv11)-indecdf(cu22,cv11));
% cc = qc;
for i = de_num
    recu22(i,:) = cu_pi(i,(seLabel(i,de_num)+1))';
    recv22(i,:) = diag(cu_pi(de_num,(reLabel(i,de_num)+1)));
    recu11(i,:) = cu_pi(i,(seLabel(i,de_num)))';
    recv11(i,:) = diag(cu_pi(de_num,(reLabel(i,de_num))));
end    
reqc = (gumbelcdf(recu22, recv22, new_theta)-gumbelcdf(recu11, recv22, new_theta)+gumbelcdf(recu11, recv11, new_theta)-gumbelcdf(recu22,recv11,new_theta));
% qc(1:de_num, 1:de_num) = reqc;
% 
reqc(reqc<=0) = 1e-16;

recc = (gumbelcdf(recu22, recv22, ctheta)-gumbelcdf(recu11, recv22, ctheta)+gumbelcdf(recu11, recv11, ctheta)-gumbelcdf(recu22,recv11,ctheta));
recc(recc<=0) = 1e-16;

theta_ratio = prod(prod(reqc./recc));

end

