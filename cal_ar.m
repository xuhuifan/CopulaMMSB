function acceptance_ratio = cal_ar( q_pi, cpi, se_Labels, re_Labels, i, ctheta, ctheta_2, dataNum, Ni, de_num)
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
% denom1 = sum(log(gumbelcdf(u22, v22, ctheta)+gumbelcdf(u11, v11, ctheta)-gumbelcdf(u21,v21, ctheta)-gumbelcdf(u12, v12, ctheta)));
qc1 = (gumbelcdf(qu22, qv22, ctheta_2)-gumbelcdf(qu11, qv22, ctheta_2)+gumbelcdf(qu11, qv11, ctheta_2)-gumbelcdf(qu22,qv11, ctheta_2));
% qc1 = (indecdf(qu22, qv22)-indecdf(qu11, qv22)+indecdf(qu11, qv11)-indecdf(qu22,qv11));
if any(de_num==i)
    qu22 = qu_pi(i,(se_Labels(i,de_num)+1))';
    qv22 = diag(qu_pi(de_num,(re_Labels(i,de_num)+1)));
    qu11 = qu_pi(i,(se_Labels(i,de_num)))';
    qv11 = diag(qu_pi(de_num,(re_Labels(i,de_num))));
    reqc1 = (gumbelcdf(qu22, qv22, ctheta)-gumbelcdf(qu11, qv22, ctheta)+gumbelcdf(qu11, qv11, ctheta)-gumbelcdf(qu22,qv11, ctheta));
    qc1(de_num) = reqc1;
end
qc1(qc1<=0) = 1e-16;

cc1 = (gumbelcdf(cu22, cv22, ctheta_2)-gumbelcdf(cu11, cv22, ctheta_2)+gumbelcdf(cu11, cv11, ctheta_2)-gumbelcdf(cu22,cv11, ctheta_2));
% cc1 = (indecdf(cu22, cv22)-indecdf(cu11, cv22)+indecdf(cu11, cv11)-indecdf(cu22,cv11));
if any(i == de_num)
    cu22 = cu_pi(i,(se_Labels(i,de_num)+1))';
    cv22 = diag(cu_pi(de_num,(re_Labels(i,de_num)+1)));
    cu11 = cu_pi(i,(se_Labels(i,de_num)))';
    cv11 = diag(cu_pi(de_num,(re_Labels(i,de_num))));
    reqc1 = (gumbelcdf(cu22, cv22, ctheta)-gumbelcdf(cu11, cv22, ctheta)+gumbelcdf(cu11, cv11, ctheta)-gumbelcdf(cu22,cv11, ctheta));
    cc1(de_num) = reqc1;
end
cc1(cc1<=0) = 1e-16;
denom1 = qc1./cc1;
clear qc1 cc1 cu22 cv22 cu11 cv11 qu22 qu11 qv22 qv11;
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
qc2 = (gumbelcdf(qu22, qv22, ctheta_2)-gumbelcdf(qu11, qv22, ctheta_2)+gumbelcdf(qu11, qv11, ctheta_2)-gumbelcdf(qu22,qv11, ctheta_2));
% qc2 = (indecdf(qu22, qv22)-indecdf(qu11, qv22)+indecdf(qu11, qv11)-indecdf(qu22,qv11));
if any(i == de_num)
    qu22 = qu_pi(i,(se_Labels(i,de_num)+1))';
    qv22 = diag(qu_pi(de_num,(re_Labels(i,de_num)+1)));
    qu11 = qu_pi(i,(se_Labels(i,de_num)))';
    qv11 = diag(qu_pi(de_num,(re_Labels(i,de_num))));
    reqc1 = (gumbelcdf(qu22, qv22, ctheta)-gumbelcdf(qu11, qv22, ctheta)+gumbelcdf(qu11, qv11, ctheta)-gumbelcdf(qu22,qv11, ctheta));
    qc2(de_num) = reqc1;
end
qc2(qc2<=0) = 1e-16;

cc2 = (gumbelcdf(cu22, cv22, ctheta_2)-gumbelcdf(cu11, cv22, ctheta_2)+gumbelcdf(cu11, cv11, ctheta_2)-gumbelcdf(cu22,cv11, ctheta_2));
% cc2 = (indecdf(cu22, cv22)-indecdf(cu11, cv22)+indecdf(cu11, cv11)-indecdf(cu22,cv11));
if any(i == de_num)
    cu22 = cu_pi(i,(se_Labels(i,de_num)+1))';
    cv22 = diag(cu_pi(de_num,(re_Labels(i,de_num)+1)));
    cu11 = cu_pi(i,(se_Labels(i,de_num)))';
    cv11 = diag(cu_pi(de_num,(re_Labels(i,de_num))));
    reqc1 = (gumbelcdf(cu22, cv22, ctheta)-gumbelcdf(cu11, cv22, ctheta)+gumbelcdf(cu11, cv11, ctheta)-gumbelcdf(cu22,cv11, ctheta));
    cc2(de_num) = reqc1;
end
cc2(cc2<=0) = 1e-16;
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

