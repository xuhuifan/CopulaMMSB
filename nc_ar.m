function acceptance_ratio = nc_ar( q_pi, cpi, se_Labels, re_Labels, i, dataNum, Ni)
%% calcluate the acceptance ratio
% q_pi: the proposal distribution for i
% cpi : the current whole distribution
% se_Labels: \{s_{ij}\}_{i,j}
% re_Labels: \{r_{ij}\}_{i,j}
% i   : the i-th data
% ctheta: copula function parameter
% dataNum : data size
% Ni:  N_{i\cdot}

cu_pi = [zeros(dataNum, 1) cumsum(cpi, 2)]; % the cumulated value \sum_{d=1}^k \pi_{id}
cu22 = cu_pi(i,(se_Labels(i,:)+1))';
cv22 = diag(cu_pi(1:dataNum,(re_Labels(i,:)+1)));
cu11 = cu_pi(i,(se_Labels(i,:)))';
cv11 = diag(cu_pi(1:dataNum,(re_Labels(i,:))));
qu_pi = cu_pi;
% fprintf('size 2 of cu_pi is %d\n', size(cu_pi, 2)); to test size(q_pi)
% fprintf('length of q_pi is %d\n', length(q_pi));
qu_pi(i,:) = [0 cumsum(q_pi)];
qu22 = qu_pi(i,(se_Labels(i,:)+1))';
qv22 = diag(qu_pi(1:dataNum,(re_Labels(i,:)+1)));
qu11 = qu_pi(i,(se_Labels(i,:)))';
qv11 = diag(qu_pi(1:dataNum,(re_Labels(i,:))));
%  denom1 = sum(log(copulacdf('Clayton', [u22 v22], ctheta)+copulacdf('Clayton', [u11 v11], ctheta)-copulacdf('Clayton', [u12 v12], ctheta)-copulacdf('Clayton', [u21 v21], ctheta)));
% denom1 = sum(log(claytoncdf(u22, v22, ctheta)+claytoncdf(u11, v11, ctheta)-claytoncdf(u21,v21, ctheta)-claytoncdf(u12, v12, ctheta)));
qc1 = (indecdf(qu22, qv22)-indecdf(qu11, qv22)+indecdf(qu11, qv11)-indecdf(qu22,qv11));
qc1(qc1==0) = 1e-16;
cc1 = (indecdf(cu22, cv22)-indecdf(cu11, cv22)+indecdf(cu11, cv11)-indecdf(cu22,cv11));
cc1(cc1==0) = 1e-16;
denom1 = qc1./cc1;
% denom1(denom1<0) = abs(denom1(denom1 < 0 ));
% denom1(isnan(denom1)) = 1;

cu22 = diag(cu_pi(1:dataNum,(se_Labels(:,i)+1)));
cv22 = cu_pi(i,(re_Labels(:,i)+1))';
cu11 = diag(cu_pi(1:dataNum,(se_Labels(:,i))));
cv11 = cu_pi(i,(re_Labels(:,i)))';
qu22 = diag(qu_pi(1:dataNum,(se_Labels(:,i)+1)));
qv22 = qu_pi(i,(re_Labels(:,i)+1))';
qu11 = diag(qu_pi(1:dataNum,(se_Labels(:,i))));
qv11 = qu_pi(i,(re_Labels(:,i)))';
% denom2 = sum(log(copulacdf('Clayton', [u22 v22], ctheta)+copulacdf('Clayton', [u11 v11], ctheta)-copulacdf('Clayton', [u12 v12], ctheta)-copulacdf('Clayton', [u21 v21], ctheta)));
qc2 = (indecdf(qu22, qv22)-indecdf(qu11, qv22)+indecdf(qu11, qv11)-indecdf(qu22,qv11));
qc2(qc2==0) = 1e-16;
cc2 = (indecdf(cu22, cv22)-indecdf(cu11, cv22)+indecdf(cu11, cv11)-indecdf(cu22,cv11));
cc2(cc2==0) = 1e-16;
denom2 = qc2./cc2;
% denom2(isnan(denom2)) = 1;
% denom2(denom2<0) = abs(denom2(denom2 < 0 ));
denom = (sum(log(denom1))+sum(log(denom2)));

a1 = sum(log(q_pi.^Ni));

a2 = sum(log(cpi(i,:).^Ni));
if isnan(a2)
    fprintf('here here\n');
end
acceptance_ratio = exp(a2-a1+denom);

end

