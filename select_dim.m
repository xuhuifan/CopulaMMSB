function [deviance, li_jps] = select_dim(dim3)
deviance = 0;

se_Labels = dim3.seLabel;
re_Labels = dim3.reLabel;
numClass = dim3.numClass;

li_jps = 0;

tau_kl = zeros(numClass, numClass);
tau1_kl = zeros(numClass, numClass);

for k = 1:numClass
    for l=1:numClass
        [x_loc, y_loc]=find((se_Labels==(k))&(re_Labels==(l)));
        tau1_kl(k,l)=sum(diag(dim3.datas(x_loc, y_loc)));
        tau_kl(k,l) =length(x_loc);
    end
end

Nikt = zeros(dim3.dataNum, numClass);
for i=1:(dim3.dataNum)
    Nikt(i, :) = N_count(se_Labels(i,:), (re_Labels(:,i))', 1:numClass);
end


for i_dev = 1:dim3.dataNum
    for j_dev = 1:dim3.dataNum
        seL = se_Labels(i_dev, j_dev);
        reL = re_Labels(i_dev, j_dev);
        
        tau_kl(seL, reL)=tau_kl(seL, reL)-1;
        tau1_kl(seL, reL)= tau1_kl(seL, reL)-dim3.datas(i_dev,j_dev);
        tau0_kl = tau_kl-tau1_kl;
        
        if (dim3.datas(i_dev,j_dev)==1)
            like_wei = (tau1_kl+dim3.lam1)./(tau_kl+dim3.lam1+dim3.lam2);
        else
            like_wei = (tau0_kl+dim3.lam2)./(tau_kl+dim3.lam1+dim3.lam2);
        end
        wei_ij = (diag(Nikt(i_dev, :))*like_wei)*diag(Nikt(j_dev, :))/((dim3.dataNum)^2);
        deviance = deviance+log(sum(sum(wei_ij)));
 
        Nikt(i_dev, seL) = Nikt(i_dev, seL)-1;
        Nikt(j_dev, reL) = Nikt(j_dev, reL)-1;
        
        if (dim3.datas(i_dev,j_dev)==1)
            like_wei = (tau1_kl(seL, reL)+dim3.lam1)/(tau_kl(seL, reL)+dim3.lam1+dim3.lam2);
        else
            like_wei = (tau0_kl(seL, reL)+dim3.lam2)/(tau_kl(seL, reL)+dim3.lam1+dim3.lam2);
        end
        
        li_jps = li_jps + log(like_wei);
         
        Nikt(i_dev, seL) = Nikt(i_dev, seL)+1;
        Nikt(j_dev, reL) = Nikt(j_dev, reL)+1;
        
        tau_kl(seL, reL)=tau_kl(seL, reL)+1;
        tau1_kl(seL, reL)= tau1_kl(seL, reL)+dim3.datas(i_dev,j_dev);
    end
end

deviance = -2*deviance;
end