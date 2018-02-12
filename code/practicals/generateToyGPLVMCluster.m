%% Genera Toy PhD. UCM-GPLVM, as in Iwata's Paper from 2016

clear all
close all
clc

K = 3;
V = 2;
Q = 3;
alpha = 0.1;
Dd = [2 2];
Nd = [500 500];
pathData = './toyData/';
addpath(genpath('E:\Clases Doctorado\Applied Bayesian Nonparametrics\GPmat-master'));
addpath(genpath('E:\Clases Doctorado\Applied Bayesian Nonparametrics\MLtoolboxes'));
% rmpath('./private');
for D = 2:2 % for several domains experiments
    
%     Md = Md(D)*ones(1,D);
%     Nd = 200*ones(1,D);
    
    for i = 1:1
        
        %         seed = 1e5*3;
%         s = RandStream('mt19937ar','Seed',1e5*i);
%         RandStream.setGlobalStream(s);
%         [S,X,Fd] = createToyGPLVMCluster(Nd,Md,D,J,K); % For the nonlinear
        %         model
        [Sclust,Y,Fd,X] = createToyMultiviewWMMCluster(Nd,Dd,V,K,Q,alpha);
%         gscatter(X{1}(:,1),X{1}(:,2),S{1})
        
%         gsave([pathData,'synth',num2str(K),'exp',num2str(i),'D',num2str(D),'J',num2str(J),'GP09_11_17.mat'],'S','W','X','Nd','Md','D','J','K');
    end
    
end

 plotMultiViewData(Y);