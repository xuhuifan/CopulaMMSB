% 
% % addpath(genpath('C:\Documents and Settings\superuser\My Documents\Dropbox\Public\copula-mmsb\lightspeed'));
% clear all;
% %% Dynamic Infinite Mixed-Membership Relational Model Sampling
% Niteration = 100000;
% samples = cell(1, Niteration);
% %% data simulation
% dataNum = 20;
% numClass = 3;
% [datas, ith_sela, ith_rela] = dataGens(dataNum);
% % E = load('senator.mat');
% % datas = E.E;
% % [dataNum, ds, tTime] = size(datas);
% %% initialization
% 
% dim3 = dim3Ini(numClass, dataNum, datas);
% dim3.datas = datas;
% 
% % whole_se = zeros(dataNum, dataNum, tTime, Niteration/2);
% % whole_re = zeros(dataNum, dataNum, tTime, Niteration/2);
% % whole_in = zeros(1, Niteration/2);
% 
% %% Gibbs sampling loop
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % deviance_numc = zeros(1, Niteration);
% % st_jp = -inf;
% % st_seL = dim3.seLabel;
% % st_reL = dim3.reLabel;
% % selec_nite = 1;
% % cu_jps = zeros(1, Niteration);
% % 
% % st_like = -inf;
% % like_seL = dim3.seLabel;
% % like_reL = dim3.reLabel;
% % selec_like = 1;
% % cu_like = zeros(1, Niteration);
% dim3.betas = ones(1, numClass);
% dim3.alphas = 2;
% 
% dim3.cpi = dirrnd(dim3.betas, dataNum);
% 
function results = FC_MMSB(datas, numClass, Niteration, de_num)
dataNum = size(datas, 1);

dim3 = dim3Ini(numClass, dataNum, datas);
dim3.datas = datas;
dim3.betas = ones(1, numClass);
dim3.alphas = 2;
dim3.cpi = dirrnd(dim3.betas, dataNum);

dim3.ctheta = gamrnd(1, 1)+1;
dim3.ctheta_2 = gamrnd(1, 1)+1;

std_deviance = zeros(1, Niteration);
std_dims = cell(1,5);
std_likeli = -inf*ones(1,5);
select_ones = zeros(1,5);
likeli_seq = zeros(1, Niteration);

tic;
%% start the for loop
for n_ite = 1:Niteration
    dim3=copula_gibbs(dim3, de_num);

    dim3.betas = beta_re(dim3.cpi, dim3.betas, dim3.numClass, dim3.dataNum, dim3.alphas);
    
    [dim3.ctheta dim3.ctheta_2] = hyper_para(dim3.dataNum, dim3.cpi, dim3.seLabel, dim3.reLabel, dim3.ctheta, dim3.ctheta_2, de_num);
%     fprintf('alpha value is %f\n', dim3.alpha);
%     fprintf('gamma value is %f\n', dim3.gamma);
%     fprintf('ctheta value is %f\n', dim3.ctheta);
    
    %     [dim3.deviance, cu_jp, cu_likes]= gibbs_dev(dim3);
%     cu_jps(n_ite) = cu_jp;
%     cu_like(n_ite) = cu_likes;
%     if cu_jp > st_jp
%         st_jp = cu_jp;
%         st_seL = dim3.seLabel;
%         st_reL = dim3.reLabel;
%         selec_nite = n_ite;
%     end
%     if cu_likes > st_like
%         st_like = cu_likes;
%         like_seL = dim3.seLabel;
%         like_reL = dim3.reLabel;
%         selec_like = n_ite;
%     end
%     
%     deviance_numc(n_ite) = dim3.deviance;
    if mod(n_ite, 100)==0
        fprintf('the iteration time is %d\n', n_ite);
    end
    
    
    [deviance, likelihoods] = select_dim(dim3);
    spe_k = ceil(5*n_ite/Niteration);
    if likelihoods > std_likeli(spe_k)
        std_likeli(spe_k) = likelihoods;
        std_dims{spe_k} = dim3;
        select_ones(spe_k) = n_ite;
    end
    likeli_seq(n_ite) = likelihoods;
    std_deviance(n_ite) = deviance;
end
toc;

 
results.std_dims = std_dims; 
results.std_deviance = std_deviance;
results.std_likeli=std_likeli;
results.select_ones = select_ones;
results.likeli_seq = likeli_seq;


