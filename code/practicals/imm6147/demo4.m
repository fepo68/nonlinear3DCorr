% Demo testing the performance of the BCD vs IRM on a dataset generated with
% community structure according to the BCD model. The analysis should
% result in the link-prediction performance as measured by AUC being worse 
% for the IRM model compared to the BCD model

clear all;
close all;

% Generate Graph without community structure
J=500;
bp=[1 1];
bn=[10 10];
gap_prior=[5 5];
alpha=5;
type='UnDirected';
[A, Z, eta, gap, A_sorted, Z_sorted, eta_sorted]=generateBCDGraph(J,alpha,bp,bn,gap_prior,type);

% Analyze the data by 5 random initialization by the IRM and BCD models
noc=10;
opts.init_sample_iter=50;
opts.nsampleiter=50;
for rep=1:5
    % Define hold-out data for cross-validation
    [W,class]=createValidationData(A,2.5,type);

    [L_IRM{rep},cpu_time_IRM{rep},Z_IRM{rep},eta_IRM{rep},sample_IRM{rep},West_IRM{rep},par_IRM{rep}]=IRMUnipartite(A,W,noc,opts);
    [AUC_IRM(rep),TPR,FPR]=calcAUC(West_IRM{rep},A);
    noc_IRM(rep)=size(Z_IRM{rep},1);

    [L_BCD{rep},cpu_time_BCD{rep},Z_BCD{rep},eta_BCD{rep},gap_BCD{rep},sample_BCD{rep},West_BCD{rep},par_BCD{rep}]=BayesianCommunityDetection(A,W,noc,opts);
    [AUC_BCD(rep),TPR,FPR]=calcAUC(West_BCD{rep},A);
    noc_BCD(rep)=size(Z_BCD{rep},1);
end

% Display the results
fprintf('Average AUC score for IRM %1.4f , Average AUC score for BCD %1.4f \n',mean(AUC_IRM),mean(AUC_BCD))
fprintf('Average number of components for IRM %6.1f ,Average number of components for BCD %6.1f \n',mean(noc_IRM),mean(noc_BCD))


%% Display a given repeat
rep=1;

figure;
J=size(A,1);
subplot(1,3,1);
mySpyPlot(A_sorted,1000/J,Z_sorted);
title('Sorted Generated Graph','FontWeight','Bold')
subplot(1,3,2);
[A_sorted_IRM,Z_sorted_IRM,eta_sorted_IRM,perm]=sortGraphUnipartite(A,sample_IRM{rep}.MAP.Z,sample_IRM{rep}.MAP.eta);
mySpyPlot(A_sorted_IRM,1000/J,Z_sorted_IRM);
title('IRM Sorted Generated Graph','FontWeight','Bold')
subplot(1,3,3);
[A_sorted_BCD,Z_sorted_BCD,eta_sorted_BCD,perm]=sortGraphUnipartite(A,sample_BCD{rep}.MAP.Z,sample_BCD{rep}.MAP.eta);
mySpyPlot(A_sorted_BCD,1000/J,Z_sorted_BCD);
title('BCD Sorted Generated Graph','FontWeight','Bold')