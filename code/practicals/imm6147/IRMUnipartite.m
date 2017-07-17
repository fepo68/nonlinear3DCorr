function [L,cpu_time,Z,eta,sample,West,par]=IRMUnipartite(A,W,noc,opts)

% Non-parametric clustering of un-directed graphs based on collapsed Gibbs
% sampling with split-merge
%
% Usage: 
%  [L,cpu_time,Z,eta,sample,West,par]=IRMUnipartite(A,W,noc,opts)
%
% Input:
%   A       cell array of I x I  sparse adjacency matrix (can be tri-angular upper matrix, notice that no self links are allowed)
%   W       cell array of I x I  sparse missing value indicator matrix (can be tri-angular upper matrix, notice that no self links are allowed)
%   noc     number of clusters
%   opts.
%           init_sample_iter initial number of sampling iterations  (i.e. burnin)
%           nsampleiter iterations to use for predictions
%           Z           noc x I matrix of intial partitioning of the nodes
%           dZstep      number of iterations between each recorded sample
%           verbose     1: display iteration results, 0: no display
%           eta0p       1x2 vector of pseudo link counts wihtin and 
%                       between communities (default: [1 1])
%           eta0n       1x2 vector of pseudo non-link counts wihtin and 
%                       between communities and (default: [1 1])
%           method      'IRM', 'IDM', 'IHW', 'IWRM', 'IWDM', 'IWHW' (default IRM)
%           sameEta     for multiple networks if sameEta=true then eta is the
%                       same across networks if sameEta=false (default) then eta is
%                       individual for each network. 
%
% Output:
%   L           Log-likelihood function at each iteration
%   cpu_time    cpu-time cost for each iteration
%   Z           Estimated clustering assigment matrix
%   eta         Estimated between group link densitites
%   sample      sampled parameters at each dZstep iterationCluster indicator matrix of size d x I where d<=noc before
%   West        Estimated link probabilities for missing links and non-links
%   par         model parameters used
%
% 
% This software is provided as is without any kind of warranty!
%
% Written by Morten Mørup.



if nargin<4
    opts=struct;
end

% Initialize variables
if ~iscell(A)
    B=A;
    clear A;
    A{1}=B;
    clear B;
end
if ~iscell(W)
    B=W;
    clear W;
    W{1}=B;
    clear B;
end
par.N=length(A);
[I,J]=size(A{1});
par.J=J;
par.type=mgetopt(opts,'type','UnDirected');
par.method=mgetopt(opts,'method','IRM');
par.sameEta=mgetopt(opts,'sameEta',false);
Iw=cell(1,par.N);
Jw=cell(1,par.N);
West=cell(1,par.N);
sumA=0;
nnzA=0;
nnzW=0;
Atot{1}=sparse(par.J,par.J);
Wtot{1}=sparse(par.J,par.J);
for n=1:par.N    
    if strcmp(par.type,'UnDirected')
        A{n}=triu(A{n});
        W{n}=triu(W{n});
    end
    W{n}=logical(W{n});
    switch par.method
        case {'IRM','IDM','IHW'}
            if par.sameEta
                 A{n}=A{n}-A{n}.*W{n}; % Remove missing links
                 Atot{1}=Atot{1}+A{n};
                 Wtot{1}=Wtot{1}+W{n};
            else
                A{n}=logical(A{n}-A{n}.*W{n}); % Remove missing links
            end
        case {'IWRM','IWDM','IWHW'}
            if par.sameEta
                A{n}=A{n}-A{n}.*W{n}; % Remove missing links
                Atot{1}=Atot{1}+A{n};
                Wtot{1}=Wtot{1}+W{n};
            else
                A{n}=A{n}-A{n}.*W{n}; % Remove missing links
            end
            
    end
    [Iw{n},Jw{n}]=find(W{n});
    West{n}=sparse(I,J);
    sumA=sumA+sum(sum(A{n}));
    nnzA=nnzA+nnz(A{n});
    nnzW=nnzW+nnz(W{n});
end

% Initialize Z
ind=ceil(noc*rand(1,I));
ind(1:noc) = 1:noc;
Z=mgetopt(opts,'Z',accumarray([ind' (1:J)'],ones(J,1),[noc J]));
noc=size(Z,1);

% Set remaining parameters
par.alpha=mgetopt(opts,'alpha',log(J));
verbose=mgetopt(opts,'verbose',1);
switch par.method
    case {'IRM','IDM','IHW'}
        par.eta0p=mgetopt(opts,'eta0p',[1 1]); % pseudo counts of links between clusters
        par.eta0n=mgetopt(opts,'eta0n',[1 1]); % pseudo counts of non-links between clusters
    case {'IWRM','IWDM','IWHW'}
        par.eta0p=mgetopt(opts,'eta0p',[1 1]); % pseudo counts of links between clusters
        par.eta0n=mgetopt(opts,'eta0n',[1 1]); % pseudo counts of total entries between clusters
end
dZstep=mgetopt(opts,'dZstep',25); % pseudo counts of links between clusters
init_sample_iter=mgetopt(opts,'init_sample_iter',10);
nsampleiter=mgetopt(opts,'nsampleiter',10);
maxiter=init_sample_iter+nsampleiter;

sample=struct;
L=zeros(1,maxiter);
cpu_time=L;
sstep=0;

iter=0;
if verbose % Display algorithm    
    disp(['Non-parametric clustering based on the Infinite ' par.method ' for ' par.type ' graphs'])
    dheader = sprintf('%12s | %12s | %12s | %12s | %12s | %12s','Iteration','logP','dlogP/|logP|','noc','time');
    dline = sprintf('-------------+--------------+--------------+--------------+--------------');
    disp(dline);
    disp(dheader);
    disp(dline);
end

Q=-inf;
Qbest=Q;
while iter<maxiter
    iter=iter+1;
    tic;
    Qold=Q;
        
    if ~par.sameEta
        % Gibbs sampling of Z in random order        
        [Z,logP_A,logP_Z]=Gibbs_sample_ZIRM(Z,A,W,randperm(J),par);            
        
        % Zplit/merge sample Z
        [Z,logP_A,logP_Z]=split_merge_sample_Z(Z,A,W,logP_A,logP_Z,par);    
        noc=size(Z,1);        
        
        % estimate mean value of eta
        eta=estimateEta(A,W,Z,par);   
    else
        % Gibbs sampling of Z in random order    
        [Z,logP_A,logP_Z]=Gibbs_sample_ZIRM(Z,Atot,Wtot,randperm(J),par);            

        % Zplit/merge sample Z
        [Z,logP_A,logP_Z]=split_merge_sample_Z(Z,Atot,Wtot,logP_A,logP_Z,par);    
        noc=size(Z,1);        
        
        % estimate mean value of eta
        eta=estimateEta(Atot,Wtot,Z,par);        
    end    
    
    Q=logP_A+logP_Z;
    dQ=Q-Qold;    
    L(iter)=Q;
    t_iter=toc;
    cpu_time(iter)=t_iter;    
    if mod(iter,dZstep)==0 
        sstep=sstep+1;
        sample.iteration(sstep)=iter;
        sample.Z{sstep}=Z;        
        sample.eta{sstep}=eta;
    end
    if Q>Qbest
        Qbest=Q;
        sample.MAP.L=Q;
        sample.MAP.iteration=iter;
        sample.MAP.Z=Z;        
        sample.MAP.eta=eta;
    end
    if rem(iter,1)==0 && verbose        
        disp(sprintf('%12.0f | %12.4e | %12.4e | %12.0f |%12.4f ',iter,Q,dQ/abs(Q),noc,t_iter));
    end
        
    % Estimate missing link probability of sample
    if iter>maxiter-nsampleiter && nnzW>0              
        if iter==maxiter-nsampleiter+1
            disp(['Initiating estimation of missing links for the last ' num2str(nsampleiter) ' iteration(s)']);   
        end
        step=10000;
        for n=1:par.N
            val=zeros(1,length(Iw{n}));
            for k=1:ceil((length(Iw{n})/step))
                 ind=(k-1)*step+1:min([k*step, length(Iw{n})]);   
                 val(ind)=sum(Z(:,Iw{n}(ind)).*(eta(:,:,n)*Z(:,Jw{n}(ind))))+eps;
            end
            West{n}=West{n}+sparse(Iw{n},Jw{n},val,I,J)/nsampleiter;                            
        end
    end
        
end

% sort Z
[val,ind]=sort(sum(Z,2),'descend');
Z=Z(ind,:);
eta=eta(ind,:,:);
eta=eta(:,ind,:);
if verbose   
  disp(sprintf('%12.0f | %12.4e | %12.4e | %12.0f |%12.4f ',iter,Q,dQ/abs(Q),noc,t_iter));
end


% -------------------------------------------------------------------------
% Parser for optional arguments
function var = mgetopt(opts, varname, default, varargin)
if isfield(opts, varname)
    var = getfield(opts, varname); 
else
    var = default;
end
for narg = 1:2:nargin-4
    cmd = varargin{narg};
    arg = varargin{narg+1};
    switch cmd
        case 'instrset',
            if ~any(strcmp(arg, var))
                fprintf(['Wrong argument %s = ''%s'' - ', ...
                    'Using default : %s = ''%s''\n'], ...
                    varname, var, varname, default);
                var = default;
            end
        otherwise,
            error('Wrong option: %s.', cmd);
    end
end
         
%---------------------------------------
function eta=estimateEta(A,W,Z,par)
        
        N=length(A);
        noc=size(Z,1);
        [n_link,n_nonlink]=calcLinkStatistics(A,W,Z,par);
        switch par.method
           case 'IRM'
                eta=n_link./(n_link+n_nonlink);
            case 'IWRM'
                eta=n_link./(n_nonlink);
            case 'IDM'                                                
                eta=zeros(noc,noc,N);
                for n=1:length(A)
                    eta(:,:,n)=diag(n_link(1:end-1,n)./(n_link(1:end-1,n)+n_nonlink(1:end-1,n)))+n_link(end,n)./(n_link(end,n)+n_nonlink(end,n))*(ones(noc)-eye(noc));
                end                        
             case 'IWDM'
                eta=zeros(noc,noc,N);
                for n=1:length(A)
                    eta(:,:,n)=diag(n_link(1:end-1,n)./n_nonlink(1:end-1,n))+n_link(end,n)./n_nonlink(end,n)*(ones(noc)-eye(noc));
                end                        
            case 'IHW'
                eta=zeros(noc,noc,N);                        
                for n=1:length(A)
                    eta(:,:,n)=sum(n_link(1:noc,n))/sum(n_link(1:noc,n)+n_nonlink(1:noc,n))*eye(noc)+n_link(end,n)/(n_link(end,n)+n_nonlink(end,n))*(ones(noc)-eye(noc));
                end
            case 'IWHW'
                eta=zeros(noc,noc,N);                        
                for n=1:length(A)
                    eta(:,:,n)=sum(n_link(1:noc,n))/sum(n_nonlink(1:noc,n))*eye(noc)+n_link(end,n)/n_nonlink(end,n)*(ones(noc)-eye(noc));
                end
        end              
        if par.N>length(A)
           eta=repmat(eta,[1 1 par.N]); 
        end        
        
        
% -------------------------------------------------------------------------  
function [Z,logP_A,logP_Z]=split_merge_sample_Z(Z,A,W,logP_A,logP_Z,par)

    %[logP_A_t,logP_Z_t]=evalLikelihood(Z,A,W,par);   
    noc=size(Z,1);
    J=size(A{1},1);    
    
    % step 1 select two observations i and j        
    ind1=ceil(J*rand);        
    ind2=ceil((J-1)*rand);
    if ind1<=ind2
       ind2=ind2+1;
    end
    clust1=find(Z(:,ind1));
    clust2=find(Z(:,ind2));

    if clust1==clust2 % Split   
        setZ=find(sum(Z([clust1 clust2],:)));    
        setZ=setdiff(setZ,[ind1,ind2]);
        n_setZ=length(setZ);
        Z_t=Z;
        Z_t(clust1,:)=0;        
        comp=[clust1 noc+1];               
        Z_t(comp(1),ind1)=1;
        Z_t(comp(2),ind2)=1;

        % Reassign by restricted gibbs sampling        
        if n_setZ>0
            for rep=1:3
                [Z_t,logP_A_t,logP_Z_t,logQ_trans,comp]=Gibbs_sample_ZIRM(Z_t,A,W,setZ(randperm(n_setZ)),par,comp);                        
            end     
        else
           logQ_trans=0;
           [logP_A_t,logP_Z_t]=evalLikelihood(Z_t,A,W,par);                 
        end

        % Calculate Metropolis-Hastings ratio
        a_split=rand<exp(logP_A_t+logP_Z_t-logP_A-logP_Z-logQ_trans);         
        if a_split
           disp(['Splitting cluster ' num2str(clust1)])
           logP_A=logP_A_t;
           logP_Z=logP_Z_t;
           Z=Z_t;
        end
    else % Merge                                     
        Z_t=Z;
        Z_t(clust1,:)=Z_t(clust1,:)+Z_t(clust2,:);
        setZ=find(Z_t(clust1,:));           
        Z_t(clust2,:)=[];        
        if clust2<clust1
            clust1_t=clust1-1;
        else 
            clust1_t=clust1;
        end
        noc_t=noc-1;

        % calculate likelihood of merged cluster       
        [logP_A_t,logP_Z_t]=evalLikelihood(Z_t,A,W,par);                

        % Zplit the merged cluster and calculate transition probabilties                
        % noc_tt=noc_t-1;
        setZ=setdiff(setZ,[ind1,ind2]);
        n_setZ=length(setZ);
        Z_tt=Z_t;
        Z_tt(clust1_t,:)=0;        
        comp=[clust1_t noc_t+1];               
        Z_tt(comp(1),ind1)=1;
        Z_tt(comp(2),ind2)=1;                
        
        % Reassign by restricted gibbs sampling
        if n_setZ>0
            for rep=1:2        
                [Z_tt,logP_A_tt,logP_Z_tt,logQ_trans,comp]=Gibbs_sample_ZIRM(Z_tt,A,W,setZ(randperm(n_setZ)),par,comp);               
            end
            Force=[1 2]*Z([clust1 clust2],:);        
            [Z_tt,logP_A_tt,logP_Z_tt,logQ_trans]=Gibbs_sample_ZIRM(Z_tt,A,W,setZ(randperm(n_setZ)),par,comp,Force);                        
        else
            logQ_trans=0;                   
        end
        a_merge=rand<exp(logP_A_t+logP_Z_t-logP_A-logP_Z+logQ_trans);                 
        
        if a_merge
          disp(['Merging cluster ' num2str(clust1) ' with cluster ' num2str(clust2)])
          logP_A=logP_A_t;
          logP_Z=logP_Z_t;
          Z=Z_t;          
        end        
    end    

%--------------------------------------------------------------------------
function [n_link,n_nonlink]=calcLinkStatistics(A,W,Z,par)

type=par.type;
method=par.method;
eta0p=par.eta0p;
eta0n=par.eta0n;
noc=size(Z,1);
N=length(A);
sumZ=sum(Z,2);
ZZT=sumZ*sumZ'-diag(sumZ);
ZZT2=sumZ.*(sumZ-1);
I=sum(sumZ);
[dummy,idx]=find(Z);
if par.sameEta
    ZZT=par.N*ZZT;
    ZZT2=par.N*ZZT2;
    ww=par.N;
else
    ww=1;
end

    switch method
        case 'IRM'
            Ap=eta0p(1)*eye(noc)+eta0p(2)*(ones(noc)-eye(noc));
            An=eta0n(1)*eye(noc)+eta0n(2)*(ones(noc)-eye(noc));       
            n_link=zeros(noc,noc,N);
            n_nonlink=zeros(noc,noc,N);
            for n=1:N
                ZAZt=Z*A{n}*Z';
                ZWZt=Z*W{n}*Z';      
                if strcmp(type,'UnDirected')   
                    n_link_tmp=ZAZt+ZAZt';    
                    n_link_tmp=n_link_tmp-0.5*diag(diag(n_link_tmp));                    
                    n_nonlink_tmp=ZZT-ZAZt-ZWZt-ZAZt'-ZWZt';
                    n_nonlink_tmp=n_nonlink_tmp-0.5*diag(diag(n_nonlink_tmp));
                else
                    n_link_tmp=ZAZt;
                    n_nonlink_tmp=ZZT-ZAZt-ZWZt;                                                                                 
                end
                n_link(:,:,n)=n_link_tmp+Ap;                           
                n_nonlink(:,:,n)=n_nonlink_tmp+An; 
            end                                        
        case 'IDM'
            n_link=zeros(noc+1,N);
            n_nonlink=zeros(noc+1,N);
            for n=1:N
                ZAZt=sum((Z*A{n}).*Z,2);
                ZWZt=sum((Z*W{n}).*Z,2);  
                if length(idx)~=par.J
                    n_link(:,n)=[ZAZt+eta0p(1); sum(sum(Z*A{n}*Z'))-sum(ZAZt)+eta0p(2)];                    
                else
                    n_link(:,n)=[ZAZt+eta0p(1); sum(sum(A{n}))-sum(ZAZt)+eta0p(2)];                    
                end
                if strcmp(type,'UnDirected')          
                    if length(idx)~=par.J
                        n_nonlink(:,n)=[ZZT2/2-ZAZt-ZWZt+eta0n(1); ww*I*(I-1)/2-sum(ZZT2)/2-(sum(sum(Z*A{n}*Z'))-sum(ZAZt))-(sum(sum(Z*W{n}*Z'))-(sum(ZWZt)))+eta0n(2)];                                        
                    else
                        n_nonlink(:,n)=[ZZT2/2-ZAZt-ZWZt+eta0n(1); ww*I*(I-1)/2-sum(ZZT2)/2-(sum(sum(A{n}))-sum(ZAZt))-(sum(sum(W{n}))-(sum(ZWZt)))+eta0n(2)];                                        
                    end
                else   
                    if length(idx)~=par.J
                        n_nonlink(:,n)=[ZZT2-ZAZt-ZWZt+eta0n(1); ww*I*(I-1)-sum(ZZT2)-(sum(sum(Z*A{n}*Z'))-sum(ZAZt))-(sum(sum(Z*W{n}*Z'))-(sum(ZWZt)))+eta0n(2)];                                       
                    else
                        n_nonlink(:,n)=[ZZT2-ZAZt-ZWZt+eta0n(1); ww*I*(I-1)-sum(ZZT2)-(sum(sum(A{n}))-sum(ZAZt))-(sum(sum(W{n}))-(sum(ZWZt)))+eta0n(2)];                                       
                    end
                end
            end
        case 'IHW'                    
            n_link=zeros(noc+1,N);
            n_nonlink=zeros(noc+1,N);
            for n=1:N
                ZAZt=sum((Z*A{n}).*Z,2);
                ZWZt=sum((Z*W{n}).*Z,2);
                if length(idx)~=par.J
                    n_link(:,n)=[ZAZt; sum(sum(Z*A{n}*Z'))-sum(ZAZt)+eta0p(2)];                    
                else
                    n_link(:,n)=[ZAZt; sum(sum(A{n}))-sum(ZAZt)+eta0p(2)];                    
                end
                if strcmp(type,'UnDirected')   
                    if length(idx)~=par.J                        
                        n_nonlink(:,n)=[ZZT2/2-ZAZt-ZWZt; ww*I*(I-1)/2-sum(ZZT2)/2-(sum(sum(Z*A{n}*Z'))-sum(ZAZt))-(sum(sum(Z*W{n}*Z'))-(sum(ZWZt)))+eta0n(2)];                                        
                    else
                        n_nonlink(:,n)=[ZZT2/2-ZAZt-ZWZt; ww*I*(I-1)/2-sum(ZZT2)/2-(sum(sum(A{n}))-sum(ZAZt))-(sum(sum(W{n}))-(sum(ZWZt)))+eta0n(2)];                                        
                    end
                else   
                    if length(idx)~=par.J
                        n_nonlink(:,n)=[ZZT2-ZAZt-ZWZt; ww*I*(I-1)-sum(ZZT2)-(sum(sum(Z*A{n}*Z'))-sum(ZAZt))-(sum(sum(Z*W{n}*Z'))-(sum(ZWZt)))+eta0n(2)];                                       
                    else
                        n_nonlink(:,n)=[ZZT2-ZAZt-ZWZt; ww*I*(I-1)-sum(ZZT2)-(sum(sum(A{n}))-sum(ZAZt))-(sum(sum(W{n}))-(sum(ZWZt)))+eta0n(2)];                                       
                    end
                end
                n_link(1,:)=n_link(1,:)+eta0p(1);
                n_nonlink(1,:)=n_nonlink(1,:)+eta0n(1);
            end
        case 'IWRM'
            n_link=zeros(noc,noc,N);
            n_nonlink=zeros(noc,noc,N);
            Ap=eta0p(1)*eye(noc)+eta0p(2)*(ones(noc)-eye(noc));
            An=eta0n(1)*eye(noc)+eta0n(2)*(ones(noc)-eye(noc));       
            for n=1:N
                ZAZt=Z*A{n}*Z';
                ZWZt=Z*W{n}*Z';  
                if strcmp(type,'UnDirected')                               
                    n_link_tmp=ZAZt+ZAZt';    
                    n_link_tmp=n_link_tmp-0.5*diag(diag(n_link_tmp));
                    n_link(:,:,n)=n_link_tmp+Ap;
                    n_nonlink_tmp=ZZT-ZWZt-ZWZt';
                    n_nonlink_tmp=n_nonlink_tmp-0.5*diag(diag(n_nonlink_tmp));
                    n_nonlink(:,:,n)=n_nonlink_tmp+An;
                else                            
                    n_link(:,:,n)=ZAZt+Ap;
                    n_nonlink(:,:,n)=ZZT-ZWZt+An;                                                          
                end
            end
        case 'IWDM'
            n_link=zeros(noc+1,N);
            n_nonlink=zeros(noc+1,N);
            for n=1:N
                ZAZt=sum((Z*A{n}).*Z,2);
                ZWZt=sum((Z*W{n}).*Z,2); 
                if length(idx)~=par.J
                    n_link(:,n)=[ZAZt+eta0p(1); sum(sum(Z*A{n}*Z'))-sum(ZAZt)+eta0p(2)];                    
                else
                    n_link(:,n)=[ZAZt+eta0p(1); sum(sum(A{n}))-sum(ZAZt)+eta0p(2)];                    
                end
                if strcmp(type,'UnDirected')                                                                 
                    if length(idx)~=par.J
                        n_nonlink(:,n)=[ZZT2/2-ZWZt+eta0n(1); ww*I*(I-1)/2-sum(ZZT2)/2-(sum(sum(Z*W{n}*Z'))-(sum(ZWZt)))+eta0n(2)];                                        
                    else
                        n_nonlink(:,n)=[ZZT2/2-ZWZt+eta0n(1); ww*I*(I-1)/2-sum(ZZT2)/2-(sum(sum(W{n}))-(sum(ZWZt)))+eta0n(2)];                                        
                    end
                else                     
                    if length(idx)~=par.J
                        n_nonlink(:,n)=[ZZT2-ZWZt+eta0n(1); ww*I*(I-1)-sum(ZZT2)-(sum(sum(Z*W{n}*Z'))-(sum(ZWZt)))+eta0n(2)];                                        
                    else
                        n_nonlink(:,n)=[ZZT2-ZWZt+eta0n(1); ww*I*(I-1)-sum(ZZT2)-(sum(sum(W{n}))-(sum(ZWZt)))+eta0n(2)];                                        
                    end
                end
            end
        case 'IWHW'
            n_link=zeros(noc+1,N);
            n_nonlink=zeros(noc+1,N);
            for n=1:N
                ZAZt=sum((Z*A{n}).*Z,2);
                ZWZt=sum((Z*W{n}).*Z,2); 
                if length(idx)~=par.J
                    n_link(:,n)=[ZAZt; sum(sum(Z*A{n}*Z'))-sum(ZAZt)+eta0p(2)];                    
                else
                    n_link(:,n)=[ZAZt; sum(sum(A{n}))-sum(ZAZt)+eta0p(2)];                    
                end
                if strcmp(type,'UnDirected')                                                                               
                    if length(idx)~=par.J
                        n_nonlink(:,n)=[ZZT2/2-ZWZt; ww*I*(I-1)/2-sum(ZZT2)/2-(sum(sum(Z*W{n}*Z'))-(sum(ZWZt)))+eta0n(2)];                                        
                    else
                        n_nonlink(:,n)=[ZZT2/2-ZWZt; ww*I*(I-1)/2-sum(ZZT2)/2-(sum(sum(W{n}))-(sum(ZWZt)))+eta0n(2)];                                        
                    end
                else                     
                    if length(idx)~=par.J
                        n_nonlink(:,n)=[ZZT2-ZWZt; ww*I*(I-1)-sum(ZZT2)-(sum(sum(Z*W{n}*Z'))-(sum(ZWZt)))+eta0n(2)];                                        
                    else
                        n_nonlink(:,n)=[ZZT2-ZWZt; ww*I*(I-1)-sum(ZZT2)-(sum(sum(W{n}))-(sum(ZWZt)))+eta0n(2)];                                        
                    end
                end                
            end
            n_link(1,:)=n_link(1,:)+eta0p(1);
            n_nonlink(1,:)=n_nonlink(1,:)+eta0n(1);
    end
    
% -------------------------------------------------------------------------
function [logP_A,logP_Z]=evalLikelihood(Z,A,W,par,n_link,n_nonlink)

    eta0p=par.eta0p;
    eta0n=par.eta0n;
    alpha=par.alpha;
    type=par.type;
    method=par.method;
    
    N=length(A);
    [I,J]=size(A{1});
    noc=size(Z,1);
    sumZ=sum(Z,2);        
    logP_A=0;    
    if nargin<6
        [n_link,n_nonlink]=calcLinkStatistics(A,W,Z,par);
    end
    
    switch method
        case 'IRM'            
            if strcmp(type,'UnDirected')            
                ii=find(repmat(triu(ones(noc)),[1,1,N])); 
                logP_A=sum(betaln(n_link(ii),n_nonlink(ii)))-N*sum([noc noc*(noc-1)/2].*betaln(eta0p,eta0n));
            else
                logP_A=sum(sum(sum(betaln(n_link,n_nonlink))))-N*sum([noc noc*(noc-1)].*betaln(eta0p,eta0n));
            end
        case 'IDM'
            logP_A=sum(sum(betaln(n_link,n_nonlink)))-N*sum([noc 1].*betaln(eta0p,eta0n));
        case 'IHW'
            link=[sum(n_link(1:end-1,:),1); n_link(end,:)];
            nolink=[sum(n_nonlink(1:end-1,:),1); n_nonlink(end,:)];
            logP_A=sum(sum(betaln(link,nolink)))-N*sum(betaln(eta0p,eta0n));
        case 'IWRM'            
            if strcmp(type,'UnDirected')            
                ii=find(repmat(triu(ones(noc)),[1,1,N])); 
                logP_A=sum(sum(clusterEvalPoisson(n_link(ii),n_nonlink(ii))))-N*sum([noc noc*(noc-1)/2].*clusterEvalPoisson(eta0p,eta0n));
            else
                logP_A=sum(sum(sum(clusterEvalPoisson(n_link,n_nonlink))))-N*sum([noc noc*(noc-1)].*clusterEvalPoisson(eta0p,eta0n));
            end
        case 'IWDM'
            logP_A=sum(sum(clusterEvalPoisson(n_link,n_nonlink)))-N*sum([noc 1].*clusterEvalPoisson(eta0p,eta0n));
        case 'IWHW'
            link=[sum(n_link(1:end-1,:),1); n_link(end,:)];
            nolink=[sum(n_nonlink(1:end-1,:),1); n_nonlink(end,:)];
            logP_A=sum(sum(clusterEvalPoisson(link,nolink)))-N*sum(clusterEvalPoisson(eta0p,eta0n));
    end    
    if nargout>1
        logP_Z=noc*log(alpha)+sum(gammaln(full(sumZ)))-gammaln(J+alpha)+gammaln(alpha);    
    end
    
% -------------------------------------------------------------------------    
function  B=clusterEvalPoisson(Np,Nn)

    B=gammaln(Np)-Np.*log(Nn);

% -------------------------------------------------------------------------
function [Z,logP_A,logP_Z,logQ_trans,comp]=Gibbs_sample_ZIRM(Z,A,W,JJ,par,comp,Force)        
    
    alpha=par.alpha;
    eta0p=par.eta0p;
    eta0n=par.eta0n;
    method=par.method;
    type=par.type;
    if nargin<7
        Force=[];
    end
    if nargin<6
        comp=[];
    end
    logQ_trans=0;

    switch method
        case {'IRM','IDM','IHW'}
            clustFun=@betaln;
        case {'IWRM','IWDM','IWHW'}
            clustFun=@clusterEvalPoisson;
    end
    
    N=length(A);
    [I,J]=size(A{1});
    eN=ones(1,N);    
    t=0;   
    sumZ=sum(Z,2);
    noc=length(sumZ);    
    q=clustFun(eta0p,eta0n);
    diag_const=q(1);
    off_const=q(2);
                
    [n_link,n_nonlink]=calcLinkStatistics(A,W,Z,par);
            
    for n=1:N
        if strcmp(type,'UnDirected')
            switch method
                case {'IRM','IDM','IHW'}                    
                    if par.sameEta
                        A{n}=A{n}+A{n}';
                    else
                        A{n}=logical(A{n}+A{n}');
                    end                    
                case {'IWRM','IWDM','IWHW'}                    
                        A{n}=A{n}+A{n}';                    
            end    
            if par.sameEta
                W{n}=W{n}+W{n}';               
            else
                W{n}=logical(W{n}+W{n}');               
            end
        else
            A2{n}=A{n}';
            W2{n}=W{n}';            
        end
    end
            
    switch method
        case {'IRM','IWRM'}
            cluster_eval=clustFun(n_link,n_nonlink);
        case {'IDM','IWDM'}            
            cluster_eval=clustFun(n_link(1:noc,:),n_nonlink(1:noc,:));              
    end
    sum_cluster_eval=zeros(1,N);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Main loop
    for k=JJ                                  
        t=t+1;
        if mod(t,5000)==0
            disp(['sampling ' num2str(t) ' out of ' num2str(J) ' nodes']);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Remove effect of s_k
        ZA1k=zeros(noc,N);
        ZW1k=zeros(noc,N);
        if ~strcmp('type','UnDirected')
            ZA2k=zeros(noc,N);
            ZW2k=zeros(noc,N);
        end
        for n=1:N
            ZA1k(:,n)=Z*A{n}(:,k);
            ZW1k(:,n)=Z*W{n}(:,k);            
            if ~strcmp(type,'UnDirected')
                ZA2k(:,n)=Z*A2{n}(:,k);                
                ZW2k(:,n)=Z*W2{n}(:,k);                                            
            end
        end        
        sumZ=sumZ-Z(:,k); 
        if par.sameEta
            switch method
                case {'IRM','IDM','IHW'}
                    nZA1k=par.N*sumZ*eN-ZA1k-ZW1k;                
                    if ~strcmp(type,'UnDirected')
                        nZA2k=par.N*sumZ*eN-ZA2k-ZW2k;
                    end
                 case {'IWRM','IWDM','IWHW'}
                    nZA1k=par.N*sumZ*eN-ZW1k;                
                    if ~strcmp(type,'UnDirected')
                        nZA2k=par.N*sumZ*eN-ZW2k;
                    end
            end
        else
            switch method
                case {'IRM','IDM','IHW'}
                    nZA1k=sumZ*eN-ZA1k-ZW1k;                
                    if ~strcmp(type,'UnDirected')
                        nZA2k=sumZ*eN-ZA2k-ZW2k;
                    end
                 case {'IWRM','IWDM','IWHW'}
                    nZA1k=sumZ*eN-ZW1k;                
                    if ~strcmp(type,'UnDirected')
                        nZA2k=sumZ*eN-ZW2k;
                    end
            end
        end
        d=find(Z(:,k));        
        % Remove link counts generated from assigment Z(:,k)
        if ~isempty(d)
             switch method
                case {'IRM','IWRM'}
                    if strcmp(type,'UnDirected')
                        n_link(:,d,:)=permute(n_link(:,d,:),[1 3 2])-ZA1k;        
                        if N==1
                            n_link(d,:)=n_link(d,:)-ZA1k';
                        else
                            n_link(d,:,:)=permute(n_link(d,:,:),[2 3 1])-ZA1k;
                        end
                        n_nonlink(:,d,:)=permute(n_nonlink(:,d,:),[1 3 2])-nZA1k;               
                        if N==1
                            n_nonlink(d,:)=n_nonlink(d,:)-nZA1k';                                      
                        else
                            n_nonlink(d,:,:)=permute(n_nonlink(d,:,:),[2 3 1])-nZA1k;                                      
                        end
                        n_link(d,d,:)=permute(n_link(d,d,:),[3 1 2])+ZA1k(d,:)';   
                        n_nonlink(d,d,:)=permute(n_nonlink(d,d,:),[3 1 2])+nZA1k(d,:)';
                    else
                        n_link(:,d,:)=permute(n_link(:,d,:),[1 3 2])-ZA1k;        
                        if N==1
                            n_link(d,:)=n_link(d,:)-ZA2k';
                        else                    
                            n_link(d,:,:)=permute(n_link(d,:,:),[2 3 1])-ZA2k;
                        end
                        n_nonlink(:,d,:)=permute(n_nonlink(:,d,:),[1 3 2])-nZA1k;               
                        if N==1
                            n_nonlink(d,:)=n_nonlink(d,:)-nZA2k';                                                                     
                        else
                            n_nonlink(d,:,:)=permute(n_nonlink(d,:,:),[2 3 1])-nZA2k;                                                                     
                        end
                    end
                case {'IDM','IWDM','IHW','IWHW'}
                    if strcmp(type,'UnDirected')
                        n_link(d,:)=n_link(d,:)-ZA1k(d,:);        
                        n_link(noc+1,:)=n_link(noc+1,:)-(sum(ZA1k,1)-ZA1k(d,:));
                        n_nonlink(d,:)=n_nonlink(d,:)-nZA1k(d,:);               
                        n_nonlink(noc+1,:)=n_nonlink(noc+1,:)-(sum(nZA1k,1)-nZA1k(d,:));
                    else
                        n_link(d,:)=n_link(d,:)-ZA1k(d,:)-ZA2k(d,:);        
                        n_link(noc+1,:)=n_link(noc+1,:)-(sum(ZA1k,1)-ZA1k(d,:))-(sum(ZA2k,1)-ZA2k(d,:));                                                                   
                        n_nonlink(d,:)=n_nonlink(d,:)-nZA1k(d,:)-nZA2k(d,:);               
                        n_nonlink(noc+1,:)=n_nonlink(noc+1,:)-(sum(nZA1k,1)-nZA1k(d,:))-(sum(nZA2k,1)-nZA2k(d,:));                                                                   
                    end                    
            end
        end
        Z(:,k)=0;               
        
        if isempty(comp) % Distinguish between restricted and non-restricted sampling
            % Non-restricted sampling
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Remove singleton cluster            
            if sumZ(d)==0 
                v=1:noc;
                v(d)=[];
                length_d=length(d);
                d=[];
                noc=noc-length_d;                               
                P=sparse(1:noc,v,ones(1,noc),noc,noc+length_d);                        
                ZA1k=P*ZA1k;
                if ~strcmp(type,'UnDirected')
                    ZA2k=P*ZA2k; 
                end
                nZA1k=P*nZA1k; 
                if ~strcmp(type,'UnDirected')
                    nZA2k=P*nZA2k; 
                end
                Z=P*Z;                        
                sumZ=sumZ(v,1);            
                switch method
                    case {'IRM','IWRM'}
                        n_link=n_link(v,v,:);            
                        n_nonlink=n_nonlink(v,v,:);            
                        cluster_eval=cluster_eval(v,v,:);            
                    case {'IDM','IWDM'}
                        n_link=n_link([v end],:);            
                        n_nonlink=n_nonlink([v end],:);            
                        cluster_eval=cluster_eval(v,:);            
                    case {'IHW','IWHW'}
                        n_link=n_link([v end],:);            
                        n_nonlink=n_nonlink([v end],:);            
                        if isempty(intersect(v,1))
                            n_link(1,:)=n_link(1,:)+eta0p(1);
                            n_nonlink(1,:)=n_nonlink(1,:)+eta0n(1);
                        end
                end
                
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Calculate probability for existing communties as well as proposal cluster            
            
            % Update cluster_eval without current node being assigned to any
            % clusters
            if ~isempty(d)
                switch method
                    case {'IRM','IWRM'}
                        cluster_eval(:,d,:)=clustFun(n_link(:,d,:),n_nonlink(:,d,:)); % removed the constant -betaln(Ap,An)))                               
                        if strcmp(type,'UnDirected')
                             cluster_eval(d,:,:)=permute(cluster_eval(:,d,:),[2 1 3]);
                        else
                            cluster_eval(d,:,:)=clustFun(n_link(d,:,:),n_nonlink(d,:,:)); % removed the constant -betaln(Ap,An)))                                 
                        end                    
                    case {'IDM','IWDM'}
                        cluster_eval(d,:)=clustFun(n_link(d,:),n_nonlink(d,:)); % removed the constant -betaln(Ap,An)))                                                       
                end                
            end
            % Evaluate total likelihood for model without the node being
            % assigned
            for n=1:N
                switch method
                    case {'IRM','IWRM'}
                        if strcmp(type,'UnDirected')
                            sum_cluster_eval(n)=sum(sum(triu(cluster_eval(:,:,n))));                        
                        else
                            sum_cluster_eval(n)=sum(sum(cluster_eval(:,:,n)));                        
                        end
                    case {'IDM','IWDM'}                        
                        sum_cluster_eval(n)=sum(cluster_eval(:,n),1);                                                                                                    
                end
            end
            
            e=ones(noc+1,1);       
            enoc=ones(1,noc);
            switch method
                case {'IRM','IWRM'}
                    if strcmp(type,'UnDirected')                
                        % Evaluate likelihood without contribution of d^th cluster
                        if N==1
                            sum_cluster_eval_d=e*sum_cluster_eval-[sum(cluster_eval,1)'; zeros(1,n)];
                        else
                            sum_cluster_eval_d=e*sum_cluster_eval-[permute(sum(cluster_eval,1),[2 3 1]); zeros(1,n)];
                        end
                        % Update likelihood when assigning node to each of the
                        % clusters
                        link=zeros(noc,noc+1,N);
                        link(:,1:noc,:)=n_link+permute(ZA1k(:,:,enoc),[1 3 2]);
                        link(:,noc+1,:)=ZA1k+eta0p(2);
                        nolink=zeros(noc,noc+1,N);
                        nolink(:,1:noc,:)=n_nonlink+permute(nZA1k(:,:,enoc),[1 3 2]);
                        nolink(:,noc+1,:)=nZA1k+eta0n(2);                
                        cluster_eval_d=clustFun(link,nolink);  
                        if N==1
                            sbeta=sum(cluster_eval_d,1)';
                        else
                            sbeta=permute(sum(cluster_eval_d,1),[2 3 1]);
                        end
                        sbeta(noc+1,:)=sbeta(noc+1,:)-noc*off_const;
                        logQ=sum(sum_cluster_eval_d+sbeta,2); % removed the constant -betaln(Ap,An)))    
                    else
                        dbeta=zeros(noc,N);
                        for n=1:N
                           dbeta(:,n)=diag(cluster_eval(:,:,n)); 
                        end
                        if N==1
                            sum_cluster_eval_d=e*sum_cluster_eval-[sum(cluster_eval+cluster_eval',1)'-dbeta; zeros(1,n)];
                        else
                            sum_cluster_eval_d=e*sum_cluster_eval-[permute(sum(cluster_eval+permute(cluster_eval,[2 1 3]),1),[2 3 1])-dbeta; zeros(1,n)];
                        end
                        
                        link1=zeros(noc,noc+1,N);
                        nolink1=zeros(noc,noc+1,N);
                        link2=zeros(noc,noc+1,N);
                        nolink2=zeros(noc,noc+1,N);                
                        for n=1:N
                            link1(:,1:noc,n)=n_link(:,:,n)+ZA1k(:,n*enoc)+diag(ZA2k(:,n));                             
                            nolink1(:,1:noc,n)=n_nonlink(:,:,n)+nZA1k(:,n*enoc)+diag(nZA2k(:,n));
                            link2(:,1:noc,n)=n_link(:,:,n)'+ZA2k(:,n*enoc)+diag(ZA1k(:,n));
                            nolink2(:,1:noc,n)=n_nonlink(:,:,n)'+nZA2k(:,n*enoc)+diag(nZA1k(:,n));
                        end
                        link1(:,noc+1,:)=ZA1k+eta0p(2);
                        nolink1(:,noc+1,:)=nZA1k+eta0n(2);                                                
                        link2(:,noc+1,:)=ZA2k+eta0p(2);                                
                        nolink2(:,noc+1,:)=nZA2k+eta0n(2);                                                                
                        cluster_eval_d=clustFun(link1,nolink1);                     
                        cluster_eval_d_2=clustFun(link2,nolink2);

                        dbeta=zeros(noc+1,N);
                        for n=1:N
                           dbeta(1:noc,n)=diag(cluster_eval_d(:,1:noc,n)); 
                        end
                        if N==1
                            sbeta=sum(cluster_eval_d+cluster_eval_d_2,1)'-dbeta;
                        else
                            sbeta=permute(sum(cluster_eval_d+cluster_eval_d_2,1),[2 3 1])-dbeta;
                        end
                        sbeta(noc+1,:)=sbeta(noc+1,:)-2*noc*off_const;
                        logQ=sum(sum_cluster_eval_d+sbeta,2); % removed the constant -betaln(Ap,An)))    
                    end 
                case {'IDM','IWDM'}
                     % Evaluate likelihood without contribution of d^th cluster                        
                        sum_cluster_eval_d=e*sum_cluster_eval-[cluster_eval; zeros(1,n)];
                        
                        % Update likelihood when assigning node to each of the
                        % clusters
                        link=zeros(noc+1,2,N);
                        nolink=link;
                        if strcmp(type,'UnDirected')
                            link(1:noc,1,:)=n_link(1:noc,:)+ZA1k;
                            link(1:noc,2,:)=n_link((noc+1)*ones(noc,1),:)+enoc'*sum(ZA1k,1)-ZA1k;
                            link(noc+1,1,:)=eta0p(1);                                                        
                            link(noc+1,2,:)=n_link(noc+1,:)+sum(ZA1k,1);
                            nolink(1:noc,1,:)=n_nonlink(1:noc,:)+nZA1k;
                            nolink(1:noc,2,:)=n_nonlink((noc+1)*ones(noc,1),:)+enoc'*sum(nZA1k,1)-nZA1k;
                            nolink(noc+1,1,:)=eta0n(1);                            
                            nolink(noc+1,2,:)=n_nonlink(noc+1,:)+sum(nZA1k,1);
                        else
                            link(1:noc,1,:)=n_link(1:noc,:)+ZA1k+ZA2k;
                            link(1:noc,2,:)=n_link((noc+1)*enoc',:)+enoc'*sum(ZA1k,1)-ZA1k+enoc'*sum(ZA2k,1)-ZA2k;
                            link(noc+1,1,:)=eta0p(1);
                            link(noc+1,2,:)=n_link((noc+1),:)+sum(ZA1k,1)+sum(ZA2k,1);
                            nolink(1:noc,1,:)=n_nonlink(1:noc,:)+nZA1k+nZA2k;
                            nolink(1:noc,2,:)=n_nonlink((noc+1)*enoc',:)+enoc'*sum(nZA1k,1)-nZA1k+enoc'*sum(nZA2k,1)-nZA2k;
                            nolink(noc+1,1,:)=eta0n(1);
                            nolink(noc+1,2,:)=n_nonlink((noc+1),:)+sum(nZA1k,1)+sum(nZA2k,1);
                        end                                                                        
                        cluster_eval_d=permute(sum(clustFun(link,nolink),2),[1 3 2]);  
                        logQ=sum(sum_cluster_eval_d+cluster_eval_d,2); % removed the constant -betaln(Ap,An)))    
                        logQ(noc+1)=logQ(noc+1)-N*diag_const;
                 case {'IHW','IWHW'}                        
                        % Update likelihood when assigning node to each of the
                        % clusters
                        link=zeros(noc+1,2,N);
                        nolink=link;                        
                        link_within=sum(n_link(1:noc,:),1); 
                        nonlink_within=sum(n_nonlink(1:noc,:),1); 
                        if strcmp(type,'UnDirected')
                            link(1:noc,1,:)=link_within(enoc',:)+ZA1k;
                            link(1:noc,2,:)=n_link((noc+1)*enoc',:)+enoc'*sum(ZA1k,1)-ZA1k;
                            link(noc+1,1,:)=link_within;                                                        
                            link(noc+1,2,:)=n_link(noc+1,:)+sum(ZA1k,1);
                            nolink(1:noc,1,:)=nonlink_within(enoc',:)+nZA1k;
                            nolink(1:noc,2,:)=n_nonlink((noc+1)*enoc',:)+enoc'*sum(nZA1k,1)-nZA1k;
                            nolink(noc+1,1,:)=nonlink_within;                            
                            nolink(noc+1,2,:)=n_nonlink(noc+1,:)+sum(nZA1k,1);
                        else
                            link(1:noc,1,:)=link_within(enoc',:)+ZA1k+ZA2k;
                            link(1:noc,2,:)=n_link((noc+1)*enoc',:)+enoc'*sum(ZA1k,1)-ZA1k+enoc'*sum(ZA2k,1)-ZA2k;
                            link(noc+1,1,:)=link_within;
                            link(noc+1,2,:)=n_link((noc+1),:)+sum(ZA1k,1)+sum(ZA2k,1);
                            nolink(1:noc,1,:)=nonlink_within(enoc',:)+nZA1k+nZA2k;
                            nolink(1:noc,2,:)=n_nonlink((noc+1)*enoc',:)+enoc'*sum(nZA1k,1)-nZA1k+enoc'*sum(nZA2k,1)-nZA2k;
                            nolink(noc+1,1,:)=nonlink_within;
                            nolink(noc+1,2,:)=n_nonlink((noc+1),:)+sum(nZA1k,1)+sum(nZA2k,1);
                        end                                                                        
                        cluster_eval_d=permute(sum(clustFun(link,nolink),2),[1 3 2]);  
                        logQ=sum(cluster_eval_d,2); % removed the constant -betaln(Ap,An)))    
                        %logQ(noc+1)=logQ(noc+1)-N*diag_const;
            end
            
                    
            % Zample from posterior     
            QQ=exp(logQ-max(logQ));
            weight=[sumZ; alpha];
            QQ=weight.*QQ;                            
            ind=find(rand<full(cumsum(QQ/sum(QQ))),1,'first');                             
            Z(ind,k)=1;   
            if ind>noc                    
                noc=noc+1;
                sumZ(noc,1)=0;
                switch method
                    case {'IRM','IWRM'}
                        n_link(:,noc,:)=eta0p(2);
                        n_link(noc,:,:)=eta0p(2);            
                        n_link(noc,noc,:)=eta0p(1);            
                        n_nonlink(:,noc,:)=eta0n(2);                     
                        n_nonlink(noc,:,:)=eta0n(2);                     
                        n_nonlink(noc,noc,:)=eta0n(1);                                 
                        cluster_eval(:,noc,:)=0;    
                        cluster_eval(noc,:,:)=0;                
                        cluster_eval_d1=permute(cluster_eval_d(:,noc,:),[1 3 2]);
                        cluster_eval_d1(noc,:)=diag_const;                                
                        ZA1k(noc,:)=0;
                        nZA1k(noc,:)=0;              
                        if ~strcmp(type,'UnDirected')
                            if N==1
                                cluster_eval_d2=cluster_eval_d_2(:,noc);
                            else
                                cluster_eval_d2=permute(cluster_eval_d_2(:,noc,:),[1 3 2]);                                                
                            end                    
                            cluster_eval_d2(noc,:)=diag_const;
                            ZA2k(noc,:)=0;
                            nZA2k(noc,:)=0;              
                            logQf=logQ(noc)+N*2*(noc-1)*off_const+N*diag_const;
                        else                
                            logQf=logQ(noc)+N*(noc-1)*off_const+N*diag_const;
                        end
                    case {'IDM','IWDM'}
                        n_link(noc+1,:)=n_link(noc,:);
                        n_link(noc,:)=eta0p(1);
                        n_nonlink(noc+1,:)=n_nonlink(noc,:);
                        n_nonlink(noc,:)=eta0n(1);                                                
                        cluster_eval(noc,:)=diag_const;                                        
                        ZA1k(noc,:)=0;                            
                        nZA1k(noc,:)=0;              
                        if ~strcmp(type,'UnDirected')                            
                            ZA2k(noc,:)=0;                            
                            nZA2k(noc,:)=0;                  
                        end
                        logQf=logQ(noc)+N*diag_const;                      
                    case {'IHW','IWHW'}
                        n_link(noc+1,:)=n_link(noc,:);
                        n_link(noc,:)=0;
                        n_nonlink(noc+1,:)=n_nonlink(noc,:);
                        n_nonlink(noc,:)=0;                                                                                                               
                        ZA1k(noc,:)=0;                            
                        nZA1k(noc,:)=0;              
                        if ~strcmp(type,'UnDirected')                            
                            ZA2k(noc,:)=0;                            
                            nZA2k(noc,:)=0;                  
                        end
                        logQf=logQ(noc);                       
                end
            else
                switch method
                    case {'IRM','IWRM'}
                        cluster_eval_d1=permute(cluster_eval_d(1:noc,ind,:),[1 3 2]);
                        if ~strcmp(type,'UnDirected')
                            cluster_eval_d2=permute(cluster_eval_d_2(1:noc,ind,:),[1 3 2]);
                        end
                end
                logQf=logQ(ind);
            end                        
        else            
            % Calculate probability for existing communties as well as proposal cluster                                                            
            switch method
                    case {'IRM','IWRM'}
                        if ~isempty(d)
                            cluster_eval(:,d,:)=clustFun(n_link(:,d,:),n_nonlink(:,d,:)); % removed the constant -betaln(Ap,An)))                               
                        end
                        if strcmp(type','UnDirected')
                            cluster_eval(d,:,:)=squeeze(cluster_eval(:,d,:));
                        else
                            cluster_eval(d,:,:)=clustFun(n_link(d,:,:),n_nonlink(d,:,:)); % removed the constant -betaln(Ap,An)))                               
                        end
                        for n=1:N
                            if strcmp(type,'UnDirected')
                                sum_cluster_eval(n)=sum(sum(triu(cluster_eval(:,:,n))));                        
                            else
                                sum_cluster_eval(n)=sum(sum(cluster_eval(:,:,n)));                        
                            end
                        end
                        e=ones(2,1);
                        if strcmp(type,'UnDirected')
                            if N==1
                                sum_cluster_eval_d=e*sum_cluster_eval-sum(cluster_eval(:,comp))';                        
                            else
                                sum_cluster_eval_d=e*sum_cluster_eval-permute(sum(cluster_eval(:,comp,:)),[2 3 1]);                        
                            end
                            link=n_link(:,comp,:)+permute(ZA1k(:,:,e),[1 3 2]);                        
                            nolink=n_nonlink(:,comp,:)+permute(nZA1k(:,:,e),[1 3 2]);            
                            cluster_eval_d1=clustFun(link,nolink);
                            if N==1
                                sbeta=sum(cluster_eval_d1,1)';            
                            else
                                sbeta=permute(sum(cluster_eval_d1,1),[2 3 1]);            
                            end
                            logQ=sum(sum_cluster_eval_d+sbeta,2); % removed the constant -betaln(Ap,An)))                                                         
                        else
                            dbeta=zeros(2,N);
                            for n=1:N
                                dbeta(:,n)=diag(cluster_eval(comp,comp,n)); 
                            end
                            if N==1
                                sum_cluster_eval_d=e*sum_cluster_eval-(sum(cluster_eval(:,comp)+cluster_eval(comp,:)')'-dbeta);                        
                            else
                                sum_cluster_eval_d=e*sum_cluster_eval-(permute(sum(cluster_eval(:,comp,:)+permute(cluster_eval(comp,:,:),[2 1 3])),[2 3 1])-dbeta);                        
                            end
                            e_noc=ones(1,2);
                            link1=n_link;
                            nolink1=n_nonlink;
                            link2=permute(n_link,[2 1 3]);
                            nolink2=permute(n_nonlink,[2 1 3]);
                            for n=1:N
                                link1(:,comp,n)=link1(:,comp,n)+ZA1k(:,n*e_noc);                                        
                                link1(comp,comp,n)=link1(comp,comp,n)+diag(ZA2k(comp,n));               
                                nolink1(:,comp,n)=nolink1(:,comp,n)+nZA1k(:,n*e_noc);                            
                                nolink1(comp,comp,n)=nolink1(comp,comp,n)+diag(nZA2k(comp));                                
                                link2(:,comp,n)=link2(:,comp,n)+ZA2k(:,n*e_noc);                                        
                                link2(comp,comp,n)=link2(comp,comp,n)+diag(ZA1k(comp,n));
                                nolink2(:,comp,n)=nolink2(:,comp,n)+nZA2k(:,n*e_noc);                            
                                nolink2(comp,comp,n)=nolink2(comp,comp,n)+diag(nZA1k(comp,n));
                            end                
                            cluster_eval_d1=clustFun(link1(:,comp,:),nolink1(:,comp,:));     
                            cluster_eval_d2=clustFun(link2(:,comp,:),nolink2(:,comp,:));     
                            dbeta=zeros(2,N);
                            for n=1:N
                               dbeta(:,n)=diag(cluster_eval_d1(comp,:,n)); 
                            end
                            if N==1
                                sbeta=sum(cluster_eval_d1+cluster_eval_d2)'-dbeta;                
                            else
                                sbeta=permute(sum(cluster_eval_d1+cluster_eval_d2,1),[2 3 1])-dbeta;                
                            end
                            logQ=sum(sum_cluster_eval_d+sbeta,2); % removed the constant -betaln(Ap,An)))                                                                             
                        end
                    case {'IDM','IWDM'}
                        link=zeros(2,N,2);
                        nolink=link;
                        if strcmp(type,'UnDirected')
                            for t=1:2                                
                                link(t,:,1)=n_link(comp(t),:)+ZA1k(comp(t),:);
                                link(t,:,2)=n_link(noc+1,:)+sum(ZA1k,1)-ZA1k(comp(t),:);
                                nolink(t,:,1)=n_nonlink(comp(t),:)+nZA1k(comp(t),:);
                                nolink(t,:,2)=n_nonlink(noc+1,:)+sum(nZA1k,1)-nZA1k(comp(t),:);
                            end
                        else                                                        
                            for t=1:2
                                link(t,:,1)=n_link(comp(t),:)+ZA1k(comp(t),:)+ZA2k(comp(t),:);
                                link(t,:,2)=n_link(noc+1,:)+sum(ZA1k,1)-ZA1k(comp(t),:)+sum(ZA2k,1)-ZA2k(comp(t),:);
                                nolink(t,:,1)=n_nonlink(comp(t),:)+nZA1k(comp(t),:)+nZA2k(comp(t),:);
                                nolink(t,:,2)=n_nonlink(noc+1,:)+sum(nZA1k,1)-nZA1k(comp(t),:)+sum(nZA2k,1)-nZA2k(comp(t),:);
                            end
                        end
                        logQ=sum(sum(clustFun(link,nolink),3),2); 
                    case {'IHW','IWHW'}
                        link=zeros(2,N,2);
                        nolink=link;                        
                        link_within=sum(n_link(1:noc,:),1);
                        nonlink_within=sum(n_nonlink(1:noc,:),1);
                        if strcmp(type,'UnDirected')
                            for t=1:2                                
                                link(t,:,1)=link_within(1,:)+ZA1k(comp(t),:);
                                link(t,:,2)=n_link(noc+1,:)+sum(ZA1k,1)-ZA1k(comp(t),:);
                                nolink(t,:,1)=nonlink_within(1,:)+nZA1k(comp(t),:);
                                nolink(t,:,2)=n_nonlink(noc+1,:)+sum(nZA1k,1)-nZA1k(comp(t),:);
                            end
                        else                                                        
                            for t=1:2
                                link(t,:,1)=link_within(1,:)+ZA1k(comp(t),:)+ZA2k(comp(t),:);
                                link(t,:,2)=n_link(noc+1,:)+sum(ZA1k,1)-ZA1k(comp(t),:)+sum(ZA2k,1)-ZA2k(comp(t),:);
                                nolink(t,:,1)=nonlink_within(1,:)+nZA1k(comp(t),:)+nZA2k(comp(t),:);
                                nolink(t,:,2)=n_nonlink(noc+1,:)+sum(nZA1k,1)-nZA1k(comp(t),:)+sum(nZA2k,1)-nZA2k(comp(t),:);
                            end
                        end
                        logQ=sum(sum(clustFun(link,nolink),3),2); 
            end
            % Zample from posterior                        
            QQ=exp(logQ-max(logQ));
            weight=sumZ(comp);
            QQ=weight.*QQ;
            QQ=QQ/sum(QQ);
            if isempty(Force)
                ind=find(rand<full(cumsum(QQ)),1,'first');
            else 
                ind=Force(k);
            end
            q_tmp=logQ-max(logQ)+log(weight);
            q_tmp=q_tmp-log(sum(exp(q_tmp)));                        
            logQ_trans=logQ_trans+q_tmp(ind);
            Z(comp(ind),k)=1;  
            switch method
                case {'IRM','IWRM'}
                    if strcmp(type,'UnDirected')
                        cluster_eval_d1=cluster_eval_d1(:,ind,:);                
                    else
                        cluster_eval_d1=permute(cluster_eval_d1(:,ind,:),[1 3 2]);
                        cluster_eval_d2=permute(cluster_eval_d2(:,ind,:),[1 3 2]);
                    end
            end            
            logQf=logQ(ind);            
            ind=comp(ind);            
        end
                        
        % Re-enter effect of new s_k        
        sumZ=sumZ+Z(:,k);
        switch method 
            case {'IRM','IWRM'}
                if strcmp(type,'UnDirected')
                    n_link(:,ind,:)=permute(n_link(:,ind,:),[1 3 2])+ZA1k;
                    if N==1
                        n_link(ind,:)=n_link(ind,:)+ZA1k';
                    else
                        n_link(ind,:,:)=permute(n_link(ind,:,:),[2 3 1])+ZA1k;
                    end
                    n_link(ind,ind,:)=permute(n_link(ind,ind,:),[3 1 2])-ZA1k(ind,:)';
                    n_nonlink(:,ind,:)=permute(n_nonlink(:,ind,:),[1 3 2])+nZA1k;        
                    if N==1
                        n_nonlink(ind,:)=n_nonlink(ind,:)+nZA1k';                
                    else
                        n_nonlink(ind,:,:)=permute(n_nonlink(ind,:,:),[2 3 1])+nZA1k;                
                    end
                    n_nonlink(ind,ind,:)=permute(n_nonlink(ind,ind,:),[3 1 2])-nZA1k(ind,:)';                
                    cluster_eval(:,ind,:)=cluster_eval_d1;
                    cluster_eval(ind,:,:)=cluster_eval_d1;
                else
                    n_link(:,ind,:)=permute(n_link(:,ind,:),[1 3 2])+ZA1k;                
                    if N==1
                        n_link(ind,:)=n_link(ind,:)+ZA2k';            
                    else
                        n_link(ind,:,:)=permute(n_link(ind,:,:),[2 3 1])+ZA2k;            
                    end
                    n_nonlink(:,ind,:)=permute(n_nonlink(:,ind,:),[1 3 2])+nZA1k;        
                    if N==1
                        n_nonlink(ind,:)=n_nonlink(ind,:)+nZA2k';                            
                    else
                        n_nonlink(ind,:,:)=permute(n_nonlink(ind,:,:),[2 3 1])+nZA2k;                            
                    end
                    cluster_eval(:,ind,:)=cluster_eval_d1;
                    if N==1
                        cluster_eval(ind,:)=cluster_eval_d2';
                    else
                        cluster_eval(ind,:,:)=cluster_eval_d2;
                    end
                end 
            case {'IDM','IWDM','IHW','IWHW'}                
                if strcmp(type,'UnDirected')
                    n_link(ind,:)=n_link(ind,:)+ZA1k(ind,:);                    
                    n_nonlink(ind,:)=n_nonlink(ind,:)+nZA1k(ind,:);                            
                    n_link(end,:)=n_link(end,:)+sum(ZA1k,1)-ZA1k(ind,:);                    
                    n_nonlink(end,:)=n_nonlink(end,:)+sum(nZA1k,1)-nZA1k(ind,:);                            
                else
                    n_link(ind,:)=n_link(ind,:)+ZA1k(ind,:)+ZA2k(ind,:);                    
                    n_nonlink(ind,:)=n_nonlink(ind,:)+nZA1k(ind,:)+nZA2k(ind,:);                          
                    n_link(end,:)=n_link(end,:)+sum(ZA1k,1)-ZA1k(ind,:)+sum(ZA2k,1)-ZA2k(ind,:);                    
                    n_nonlink(end,:)=n_nonlink(end,:)+sum(nZA1k,1)-nZA1k(ind,:)+sum(nZA2k,1)-nZA2k(ind,:);                            
                end
                if strcmp(method,'IDM') || strcmp(method,'IWDM')
                    cluster_eval(ind,:)=clustFun(n_link(ind,:),n_nonlink(ind,:));                     
                end                
        end                   
        % Remove empty clusters        
        if ~all(sumZ)
            d=find(sumZ==0);
            ind_d=find(d<comp);
            comp(ind_d)=comp(ind_d)-1;
            v=1:noc;
            v(d)=[];
            noc=noc-1;                               
            P=sparse(1:noc,v,ones(1,noc),noc,noc+1);                        
            Z=P*Z;                        
            sumZ=sumZ(v,1);    
            switch method 
                case {'IRM','IWRM'}
                    n_link=n_link(v,v,:);
                    n_nonlink=n_nonlink(v,v,:);
                    cluster_eval=cluster_eval(v,v,:);            
                case {'IDM','IWDM'}
                    n_link=n_link([v end],:);
                    n_nonlink=n_nonlink([v end],:);
                    cluster_eval=cluster_eval(v,:);            
                case {'IHW','IWHW'}
                    n_link=n_link([v end],:);
                    n_nonlink=n_nonlink([v end],:);
                    if isempty(intersect(v,1))
                        n_link(1,:)=n_link(1,:)+eta0p(1);
                        n_nonlink(1,:)=n_link(1,:)+eta0n(1);
                    end
            end            
        end           
        
    end              
    noc=length(sumZ);
    logP_Z=noc*log(alpha)+sum(gammaln(full(sumZ)))-gammaln(J+alpha)+gammaln(alpha);
    switch method
        case {'IRM','IWRM'}
            if strcmp(type,'UnDirected')
                logP_A=logQf-N*sum([noc noc*(noc-1)/2].*[diag_const off_const]);                           
            else
                logP_A=logQf-N*sum([noc noc*(noc-1)].*[diag_const off_const]);                                           
            end                
        case {'IDM','IWDM'}                     
            logP_A=sum(sum(clustFun(n_link,n_nonlink)))-N*sum([noc 1].*[diag_const off_const]);                           
        case {'IHW','IWHW'}
            logP_A=logQf-N*sum([diag_const off_const]);                           
    end
    
    %if strcmp(type,'UnDirected')
    %    for n=1:length(A)
    %        A{n}=triu(A{n});
    %        W{n}=triu(W{n});
    %    end
    %end        
    %[n_link_t,n_nonlink_t]=calcLinkStatistics(A,W,Z,par);   
    %if round(sum(n_link(:)-n_link_t(:)))~=0 || round(sum(n_nonlink(:)-n_nonlink_t(:)))~=0
    %    sum(n_link(:)-n_link_t(:))
    %    sum(n_nonlink(:)-n_nonlink_t(:))
    %    keyboard;
    %end