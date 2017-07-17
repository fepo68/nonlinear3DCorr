% This scripts analyses the synthetic example given in the paper with 5 small and 3 large communities
load DemoData2;

% Create Validation Data
pct_missing=2.5;
type='UnDirected';
[W,class]=createValidationData(A,pct_missing,type);

% Run the BCD algorithm
opts.init_sample_iter=25; % Use 25 burn in iterations
opts.nsampleiter=25;      % Sample for prediction of link for 25 iterations
noc=5;
[L,cpu_time,Z,eta,gap,sample,West,predL]=BayesianCommunityDetection(A,W,noc,opts);

% Plot the results
plotSyntheticResults(A,West,Z_true,sample.MAP.Z);

%% Run the IRM model
[L_IRM,cpu_time_IRM,Z_IRM,eta_IRM,sample_IRM,West_IRM,par_IRM]=IRMUnipartite(A,W,noc,opts);

% Plot the results
plotSyntheticResults(A,West_IRM,Z_true,sample_IRM.MAP.Z);

