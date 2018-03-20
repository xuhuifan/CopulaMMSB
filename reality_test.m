% reality mining experiment
load('reality_mit.mat');

datas = network.lab;

datas(datas <= 1) = 0;
datas(isnan(datas)) = 0;
datas(datas > 1) = 1;
%datas = network.friends;
%datas(isnan(datas)) = 0;

numClass = 4;
Niteration = 600;

 results = FNC_MMSB(datas, numClass, Niteration);

% results = FC_MMSB(datas, numClass, Niteration);

 % wiki_talk test
 clear;
 load('dophin.mat');
 numClass = 4;
 Niteration = 60000;
 results = FNC_MMSB(rela_mat, numClass, Niteration);
 save dolphis_FC_2.mat results;
 
 datas = dataGens_2(50, numClass);
 a = FNC_MMSB(datas, numClass, Niteration);
 
%   % synthetic data testing 2_11
%  clear;
%  dataNum = 50;
%  numClass = 4;
%  datas = dataGens_2(dataNum, numClass);
% reality mining experiment
load('reality_mit.mat');

datas = network.lab;

datas(datas <= 1) = 0;
datas(isnan(datas)) = 0;
datas(datas > 1) = 1;

 Niteration = 40000;
 load('de_num.mat');
% de_num = [1:2 5:20];
numClass = 4;
 a = FC_MMSB(datas, numClass, Niteration, de_num);
 
 save rea_FC_1.mat a;
 %%%%%
 clear;
 
 load('syn.mat');

 Niteration = 1000;
 de_num =  1:20;
 result = FNC_MMSB(datas, numClass, Niteration);
 
 
 
 
 