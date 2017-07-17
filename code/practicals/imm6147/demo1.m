% This scripts generates a graph according to the Bayesian Community Detection (BCD) model and infer the
% parameters of the model using BCD and IRM

J=250; % Number of nodes
alpha=5; % Parameter to Z~CRP(alpha)
bp=[1 1];
bn=[1 1]; %             eta_ll~Beta(bp(1),bn(1)), eta_lm~trucatedBeta(bp(2),bn(2))
gap_prior=[1 1]; %      gap~Beta(gap_prior(1),gap_prior(1));
type='UnDirected';

[A,Z_true,eta_true,gap]=generateBCDGraph(J,alpha,bp,bn,gap_prior,type);

% Create Validation Data
pct_missing=2.5;
[W,class]=createValidationData(A,pct_missing,type);

% Run the BCD algorithm
opts.init_sample_iter=25; % Use 25 burn in iterations
opts.nsampleiter=25;      % Sample for prediction of link for 25 iterations
noc=5;
[L,cpu_time,Z,eta,gap,sample,West,predL]=BayesianCommunityDetection(A,W,noc,opts);
% Plot the results
plotSyntheticResults(A,West,Z_true,sample.MAP.Z);

%% Run the IRM model
[L_IRM,cpu_time_IRM,Z_IRM,eta_IRM,sample_IRM,West_IRM, par]=IRMUnipartite(A,W,noc,opts);
% Plot the results
plotSyntheticResults(A,West,Z_true,sample_IRM.MAP.Z);

