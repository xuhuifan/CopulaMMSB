function dim3 = copula_gibbs(dim3, de_num)
%% using block-slice sampling scheme to do the sampling
se_Labels = dim3.seLabel;  % data's sender's label
re_Labels = dim3.reLabel;  % data's receiver's label
dataNum = dim3.dataNum;    % data number
% indexLabel = dim3.indexLabel; % cluster label index
ctheta = dim3.ctheta;  % copula function's parameter
ctheta_2 = dim3.ctheta_2;
cpi = dim3.cpi; % current \{\pi_i\}_{i=1}^n

numClass = dim3.numClass;

tau_kl = dim3.tau_kl;
tau1_kl = dim3.tau1_kl;
%  fprintf('sum of tau_kl is %d\n', sum(sum(tau_kl)));  %%%%% one problem here

Nikt = zeros(dataNum, numClass); % pre_define N_{ik}

for i_sli=randperm(dataNum) % N_{ik}'s value
    %     Nikt(i_sli, 1:numClass) = N_count(se_Labels(i_sli,:), (re_Labels(:,i_sli))', indexLabel);
    Nikt(i_sli, :) = N_count(se_Labels(i_sli,:), (re_Labels(:,i_sli))', 1:numClass);
end
% fprintf('length of betas is %d\n', length(betas));
% fprintf('numClass is %d\n', numClass);

%% sample \{\pi_i\}_{i=1}^n
% q_pi=dirrnd((repmat(dim3.betas, dataNum, 1)+Nikt), dataNum);    % the proposal distribution \{q_i\}_{i=1}^n = \{pi_i^*\}_{i=1}^n

ar = zeros(dataNum, 1);
br = zeros(dataNum, 1);
for i = randperm(dataNum)
    q_pi = dirrnd((dim3.alphas*dim3.betas+Nikt(i,:)), 1);
    acceptance_ratio = cal_ar(q_pi, cpi, se_Labels, re_Labels, i, ctheta, ctheta_2, dataNum, Nikt(i,:), de_num);
    ar(i) = acceptance_ratio;
    if rand < acceptance_ratio
        cpi(i,:) = q_pi;  % update the current \pi_i
        br(i) = 1;
    end
    if isnan(acceptance_ratio)
        fprintf('NaN again! \n');
    end
end
% ar(ar>1)=1;
%  fprintf('average acceptance ratio is %f\n', mean(ar));


%% sample \{s_{ij}, r_{ij}\}_{i,j}
cu_pi = [zeros(dataNum, 1) cumsum(cpi, 2)]; % the cumulated value \sum_{d=1}^k \pi_{id}, add zeros(dataNum, 1) because we need to consider the 0 case in copula function
for i=randperm(dim3.dataNum)
    for j=randperm(dim3.dataNum)
        % here may lie some problems
        
        %         uu = repmat(cu_pi(i,:)', 1, numClass+2);  % \{u_i^k\}_{k=0}^{K+1} in Eq. (6)
        %         vv = repmat(cu_pi(j,:), numClass+2, 1);   % \{v_j^k\}_{k=0}^{K+1} in Eq. (6)
        % %% use another to calculate copula cdf value
        % %         cprob = copulacdf('Clayton', [uu(:) vv(:)], ctheta);  % calculate all the copula function value (be one long (K+2)*(K+2) column)
        % %         cprob = reshape(cprob, (numClass+2), (numClass+2));    % reshape the column to the matrix case
        %
        %         cprob = gumbelcdf(uu, vv, ctheta);
        %         % p_{ij}^{kl} value in Eq. (6) \forall k,l\in\{1,...,K+1\}
        %         prob_can = cprob(2:(numClass+2),2:(numClass+2))-cprob(1:(numClass+1),2:(numClass+2))+cprob(1:(numClass+1),1:(numClass+1))-cprob(2:(numClass+2),1:(numClass+1));
        
        u22 = repmat(cu_pi(i,2:(numClass+1))', 1, numClass);   %use another way to compute the likelihood p(s_{ij}=k, r_{ij}=l)
        v22 = repmat(cu_pi(j,2:(numClass+1)), numClass, 1);
        u11 = repmat(cu_pi(i,1:(numClass))', 1, numClass);
        v11 = repmat(cu_pi(j,1:(numClass)), numClass, 1);
        if (any(i== de_num)&any(j== de_num))
            prob_can = gumbelcdf(u22, v22, ctheta)-gumbelcdf(u11, v22, ctheta)+gumbelcdf(u11, v11, ctheta)-gumbelcdf(u22,v11, ctheta);
        else
            prob_can = gumbelcdf(u22, v22, ctheta_2)-gumbelcdf(u11, v22, ctheta_2)+gumbelcdf(u11, v11, ctheta_2)-gumbelcdf(u22,v11, ctheta_2);
        end
        prob_can(prob_can < 0) = 0;
        %         if sum(sum(prob_can<0))
        %             prob_can
        %             fprintf('here here\n');
        %         end
        
        %         prob_can = copulacdf('Gaussian', [u22(:) v22(:)], ctheta)+copulacdf('Gaussian', [u11(:) v11(:)], ctheta)-copulacdf('Gaussian', [u12(:) v12(:)], ctheta)-copulacdf('Gaussian', [u21(:) v21(:)], ctheta);
        %         prob_can = prob_can/sum(prob_can);
        %         prob_can = reshape(prob_can, (numClass+1), (numClass+1));
        
        %         se_la = find(indexLabel == se_Labels(i,j));
        %         re_la = find(indexLabel == re_Labels(i,j));
        se_la = se_Labels(i,j);
        re_la = re_Labels(i,j);
        % edge(i,j,t)'s likelihood calculation
        tau_kl(se_la, re_la)=tau_kl(se_la, re_la)-1;
        tau1_kl(se_la, re_la)= tau1_kl(se_la, re_la)-dim3.datas(i,j);
        tau0_kl = tau_kl-tau1_kl;
        
        % calculating the likehood value
        if (dim3.datas(i,j)==1)
            like_wei = (tau1_kl+dim3.lam1)./(tau_kl+dim3.lam1+dim3.lam2);  % change the denominator
        else
            like_wei = (tau0_kl+dim3.lam2)./(tau_kl+dim3.lam1+dim3.lam2);
        end
        
        p_weight = prob_can.*like_wei;
        %             end
        
        p_weight = reshape(p_weight, 1, []);
        
        % sampling
        ath_value = 1+sum((rand*sum(p_weight)) > cumsum(p_weight));
        ath_col = ceil(ath_value/numClass);
        ath_row = ath_value - (ath_col-1)*numClass;
        
        se_Labels(i,j) = (ath_row);
        re_Labels(i,j) = (ath_col);
        
        tau_kl(ath_row, ath_col)=tau_kl(ath_row, ath_col)+1;
        tau1_kl(ath_row, ath_col)=tau1_kl(ath_row, ath_col)+dim3.datas(i,j);
        
    end
end

dim3.seLabel = se_Labels;
dim3.reLabel = re_Labels;
dim3.tau_kl = tau_kl;
dim3.tau1_kl = tau1_kl;
% dim3.indexLabel = 1:numClass;
dim3.cpi = cpi;

end