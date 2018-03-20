function [datas ith_sela ith_rela] = dataGens(dataNum)

% datas = binornd(ones(datanum, datanum, ttime),pvalue*ones(datanum, datanum, ttime));
simi_mat = [0.95 0.05 0 0;0.05 0.95 0.05 0; 0.05 0 0.95 0;0 0.05 0 0.95];
ith_sela = zeros(dataNum, dataNum);
ith_rela = zeros(dataNum, dataNum);

mm_value = [0.8 0.2 0 0; 0.1 0.7 0.1 0.1; 0.1 0.05 0.8 0.05; 0.1 0.1 0.1 0.7];
num_cl = size(mm_value, 1);
de_num = [ceil(dataNum/num_cl*2)-5 ceil(dataNum/num_cl) ceil(dataNum/num_cl)-5 dataNum-2*ceil(dataNum/num_cl)-ceil(dataNum/num_cl*2)+10];
de_nums = cumsum(de_num);
k_clu = zeros(1,dataNum);
for i=1:num_cl
    if i > 1
        k_clu((de_nums(i-1)+1):de_nums(i)) = i;
    else
        k_clu(1:de_nums(i)) = i;
    end
end

for i=1:dataNum
    for j=1:dataNum
        ith_sela(i, j) = 1+sum(rand > cumsum(mm_value(k_clu(i), :)));
        ith_rela(i, j) = 1+sum(rand > cumsum(mm_value(k_clu(j), :)));
    end
end

    for ith = 1:dataNum
        for jth = 1:dataNum
            datas(ith,jth) = binornd(1, simi_mat(ith_sela(ith,jth), ith_rela(ith,jth)));
        end
    end