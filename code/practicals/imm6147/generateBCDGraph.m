function [A,Z,eta, gap, A_sorted,Z_sorted,eta_sorted]=generateBCDGraph(J,alpha,bp,bn,gap_prior,type)
% Script to generate binary unipartite graphs according to the Bayesian Community Detection
% (BCD) generative model
%
% Usage:
%   [A,Z,eta,gap]=generateBCDGraph(J,alpha,bp,bn,gap_par)
% 
% Input:
%   J           size of graph 
%   alpha       parameter for the Chinese Restaurant Process (CRP) 
%   bp          1x2 vector (default [10 1]) where 
%               bp(1): indicate within community prior link count
%               bp(2): indicate between community prior link count
%   bn          1x2 vector (default [1 10]) where 
%               bn(1): indicate within community prior non-link count
%               bn(2): indicate between community prior non-link count
%   gap_prior   prior distribution parameters for the gap parameter (default [1 1])
%   type        'UnDirected' (default) or 'Directed'
%
% Output:
%   A           Generated graph that has not been sorted, i.e. A=Bernoulli(Z'*eta*Z)
%   Z           noc x J generated assignment matrix, i.e. Z ~ Discrete(mu)
%   eta         noc x noc generated group relations, i.e. eta_{lm} ~ Beta(Bp,Bn)
%   gap         gap parameter
%   A_sorted    sorted J x J adjacency matrix of the generated graph,
%   Z_sorted    sorted noc x J generated assignment matrix
%   eta_sorted  sorted noc x noc generated group relations
%   perm        permutation matrix, i.e. A_sorted=A(perm,perm)
%
% Written by Morten Mørup

if nargin<6
    type='UnDirected';
end
if nargin<5
    gap_prior=[1 10];
end
if nargin<4
    bn=[1 10];    
end
if nargin<3
    bp=[10 1];
end
if nargin<2
    alpha=1;
end


%%%%%%%%%%%%%%%%%%
% Generate graph %
%%%%%%%%%%%%%%%%%%
% We Draw the columns z_i of the assignemnt matrix Z from the Chinese
% Restaurant Process
Z=zeros(1,J);
Z(1,1)=1;
sumZ=sum(Z,2);
for i=2:J
   p=[sumZ alpha]./(alpha+i-1);
   pp=cumsum(p);
   ind=find(rand<pp,1,'first');
   Z(ind,i)=1;      
   if ind>length(sumZ)
       sumZ=[sumZ 1];
   else
       sumZ(ind)=sumZ(ind)+1;
   end
end
noc=size(Z,1);

% Parameter for the beta prior imposed on eta
% We will assume the links drawn within and between communities are drawn
% from the same distribution specified by bp(1), bn(1) and bp(2), bn(2) respectively

% We next draw cluster relations eta from the Beta distribution specified
% by a and b
% Draw within community densities
eta_diag=betarnd(bp(1)*ones(noc,1),bn(1)*ones(noc,1));
gap=betarnd(gap_prior(1),gap_prior(2));

% Draw between community densities
eta=zeros(noc);
for i=1:noc
    for ii=i+1:noc
        p=betainc(gap*min([eta_diag(i), eta_diag(ii)]),bp(2),bn(2));
        p=p.*rand;
        eta(i,ii) = betainv(p,bp(2),bn(2));
    end
end
if strcmp(type,'UnDirected')
    eta=eta+eta'+diag(eta_diag);
end

% We finally draw links A from the Bernoulli likelihood
% We are currently interested in drawing an undirected graph thus we draw half of the links from
% the upper and half of the links from the lower triangular part of the
% adjacency matrix and add these parts together.
A = (Z'*eta*Z)>rand(J);
if strcmp(type,'UnDirected')
    A=triu(A,1);
    A=A+A';
    A=sparse(A);
end

[A_sorted,Z_sorted,eta_sorted]=sortGraphUnipartite(A,Z,eta,ones(1,J));

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display generated graph %
%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
subplot(2,2,1);
mySpyPlot(A,[],1000/J);
title('Generated Graph','FontWeight','Bold')
subplot(2,2,2);
[A_sorted,Z_sorted,eta_sorted,perm]=sortGraphUnipartite(A,Z,eta);
mySpyPlot(A_sorted,1000/J,Z_sorted);
title('Sorted Generated Graph','FontWeight','Bold')
subplot(2,2,3);
imagesc(Z_sorted); colormap(1-gray); axis off; title('Generated Sorted Clustering Assigment Matrix Z','FontWeight','Bold')
subplot(2,2,4);
imagesc(eta_sorted); colormap(1-gray); axis equal; axis tight; title('Generated Between Cluster Relations \eta','FontWeight','Bold')



