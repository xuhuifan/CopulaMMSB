function [datas ith_sela ith_rela] = dataGens_2(dataNum, numClass)

% datas = binornd(ones(datanum, datanum, ttime),pvalue*ones(datanum, datanum, ttime));
simi_mat = [0.95 0.05 0 0;0.05 0.95 0.05 0; 0.05 0 0.95 0;0 0.05 0 0.95];
ith_sela = zeros(dataNum, dataNum);
ith_rela = zeros(dataNum, dataNum);
datas = zeros(dataNum, dataNum);

mm_value = [0.9 0.1 0 0; 0 0.9 0.1 0; 0.1 0.05 0.85 0; 0.1 0.05 0.05 0.8];
num_cl = size(mm_value, 1);
de_num = [ceil(dataNum/num_cl*2)-5 ceil(dataNum/num_cl) ceil(dataNum/num_cl)-5 dataNum-2*ceil(dataNum/num_cl)-ceil(dataNum/num_cl*2)+10];
true_value = [repmat(mm_value(1,:), de_num(1), 1);repmat(mm_value(2,:), de_num(2), 1);repmat(mm_value(3,:), de_num(3), 1);repmat(mm_value(4,:), de_num(4), 1)];

cumulative_tv = [zeros(dataNum, 1) cumsum(true_value, 2)];

true_theta = 3.5;

for gen_i = 1:dataNum
    for gen_j = 1:dataNum
        u22 = repmat(cumulative_tv(gen_i,2:(numClass+1))', 1, numClass);
        v22 = repmat(cumulative_tv(gen_j,2:(numClass+1)), numClass, 1);
        u11 = repmat(cumulative_tv(gen_i,1:(numClass))', 1, numClass);
        v11 = repmat(cumulative_tv(gen_j,1:(numClass)), numClass, 1);
        if ((gen_i<= de_num(1))&(gen_i<= de_num(1)))
            p_weight = gumbelcdf(u22, v22, true_theta)-gumbelcdf(u11, v22, true_theta)+gumbelcdf(u11, v11, true_theta)-gumbelcdf(u22,v11, true_theta);
        else
            p_weight = indecdf(u22, v22)-indecdf(u11, v22)+indecdf(u11, v11)-indecdf(u22,v11);
        end
        p_weight(p_weight < 0) = 0;
        
        p_weight = reshape(p_weight, 1, []);
        
        % sampling
        ath_value = 1+sum((rand*sum(p_weight)) > cumsum(p_weight));
        ath_col = ceil(ath_value/numClass);
        ath_row = ath_value - (ath_col-1)*numClass;
        
        ith_sela(gen_i,gen_j) = (ath_row);
        ith_rela(gen_i,gen_j) = (ath_col);
        
        datas(gen_i, gen_j) = (rand > simi_mat(ath_row, ath_col));
    end
end


for ith = 1:dataNum
    for jth = 1:dataNum
        datas(ith,jth) = binornd(1, simi_mat(ith_sela(ith,jth), ith_rela(ith,jth)));
    end
end