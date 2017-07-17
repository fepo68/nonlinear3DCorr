function [A, Z, eta, A_sorted, Z_sorted, eta_sorted, perm]=generateIRMGraph(J,alpha,bp,bn,type)
% Script to generate binary unipartite graphs according to the Infinite Relational Model
%
% Usage:
%   [A,Z,eta,A_sorted,Z_sorted,eta_sorted,perm]=generateIRMGraph(J,alpha,bp,bn,type)
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
%   type        'UnDirected' (default) or 'Directed'
%
% Output:
%   A           Generated graph that has not been sorted, i.e. A=Bernoulli(Z'*eta*Z)
%   Z           noc x J generated assignment matrix, i.e. Z ~ Discrete(mu)
%   eta         noc x noc generated group relations, i.e. eta_{lm} ~ Beta(Bp,Bn)
%   A_sorted    sorted J x J adjacency matrix of the generated graph,
%   Z_sorted    sorted noc x J generated assignment matrix
%   eta_sorted  sorted noc x noc generated group relations
%   perm        permutation matrix, i.e. A_sorted=A(perm,perm)
%
% Written by Morten Mørup

if nargin<5
    type='UnDirected';
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

% Z ~ CRP(alpha)
% We Draw the columns z_i of the assignemnt matrix Z from the Chinese
% Restaurant Process
Z=zeros(1,J);
Z(1,1)=1;       % assign first customer to first table
sumZ=sum(Z,2);  % number of customers currently seated at the various tables
for i=2:J % iterate over remaining customers   
   % Draw an existing or new table from the probability distribution p
   p=[sumZ alpha]./(alpha+i-1);
   pp=cumsum(p);
   ind=find(rand<pp,1,'first');
   
   % Place customer at existing table given by ind     
   Z(ind,i)=1;    
   % Update the sum of customers at each table
   if ind>length(sumZ) 
       sumZ=[sumZ 1];
   else
       sumZ(ind)=sumZ(ind)+1; 
   end
end
noc=size(Z,1); % Number of generated clusters, i.e. non-empty tables

% eta_lm ~ beta(bp,bn)
% (We will assume the links drawn within and between communities are drawn
% from the same distribution specified by bp(1), bn(1) and bp(2), bn(2)
% respectively)
eta = betarnd(bp(1)*eye(noc)+bp(2)*(ones(noc)-eye(noc)),bn(1)*eye(noc)+bn(2)*(ones(noc)-eye(noc)));
if strcmp(type,'UnDirected') % Force eta=eta' for UnDirected graphs
     eta=triu(eta,1)+triu(eta,1)'+diag(diag(eta));
end

% A_ij ~ Bernoulli(z_i'*eta*z_j)
A = (Z'*eta*Z)>rand(J);

% Remove self-links
A=A-diag(diag(A));
if strcmp(type,'UnDirected')
    A=triu(A,1);
    A=A+A';
end
A=sparse(A);
    

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


