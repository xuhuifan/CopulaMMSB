function dim3 = dim3Ini(numClass, dataNum, datas)



dim3.dataNum = dataNum;

re_Labels = zeros(dataNum, dataNum);
se_Labels = zeros(dataNum, dataNum);

    for i=1:dataNum
        re_value = mnrnd(1, 1/numClass*ones(1, numClass), dataNum);
        se_value = mnrnd(1, 1/numClass*ones(1, numClass), dataNum);        
        for j=1:dataNum
            re_Labels(i,j) = find(re_value(j,:)==1);
            se_Labels(i,j) = find(se_value(j,:)==1);
        end
    end

    tau_kl = zeros(numClass, numClass);  % pre_define the matrix n_{k,l}
tau1_kl = zeros(numClass, numClass); % pre_define the matrix n_{k,l}^1

for k = 1:numClass    % this is to calculate the matrix n_{k,l}'s value
    for l=1:numClass
        % [x_loc, y_loc]=find((se_Labels==indexLabel(k))&(re_Labels==indexLabel(l)));
        xy_loc=((se_Labels==(k))&(re_Labels==(l)));
        tau1_kl(k,l)=sum(row_sum(datas(xy_loc)));
        tau_kl(k,l) =sum(row_sum(xy_loc));
    end
end

dim3.lam1 = sum(sum(tau1_kl))/sum(sum(tau_kl));
dim3.lam2 = 1;

dim3.numClass = numClass;

dim3.reLabel = re_Labels;
dim3.seLabel = se_Labels;

dim3.tau_kl = tau_kl;
dim3.tau1_kl = tau1_kl;
