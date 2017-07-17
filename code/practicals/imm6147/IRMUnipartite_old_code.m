function [L,cpu_time,Z,eta,sample,West,predL]=IRMUnipartite(A,W,noc,opts)

% Non-parametric clustering of un-directed graphs based on collapsed Gibbs sampling with split
% merge
%
% Usage: 
%  [L,cpu_time,Z,eta,sample,West]=IRMmisDP(A,W,noc,opts)
%
% Input:
%   A       cell array of I x I  sparse adjacency matrix (can be triuangular upper matrix)
%   W       cell array of I x I  sparse missing value indicator matrix (can be triuangular upper matrix)
%   noc     number of clusters
%   opts.
%           maxiter     maximum number of iterations
%           Z           I x noc matrix of community identities
%           dZstep      number of iteratiosn between each recorded sample
%           verbose     1: display iteration results, 0: no display
%           eta0p       1x2 vector of pseudo link counts wihtin and 
%                       between communities (default: [5 1])
%           eta0n       1x2 vector of pseudo non-link counts wihtin and 
%                       between communities and (default: [1 5]) 
%           method      'IRM', 'IDM', 'IHM', 'IWRM', 'IWDM', IWHM (default IRM)
%
% Output:
%   L           Log-likelihood function at each iteration
%   cpu_time    cpu-time cost for each iteration
%   Z           Estimated clustering assigment matrix
%   sample      sampled parameters at each dZstep iterationCluster indicator matrix of size d x I where d<=noc before
%   eta         Estimated between community link densitites
%   West        Estimated link probabilities for missing links and non-links
%   predL       predictive Log-likelihood
 
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
N=length(A);
[I,J]=size(A{1});
type=mgetopt(opts,'type','UnDirected');
method=mgetopt(opts,'method','IRM');
Iw=cell(1,N);
Jw=cell(1,N);
West=cell(1,N);
predL=zeros(1,N);
sumA=0;
nnzA=0;
for n=1:N    
    if strcmp(type,'UnDirected')
        A{n}=triu(A{n});
        W{n}=triu(W{n});
    end
    W{n}=logical(W{n});
    switch method
        case {'IRM','IDM','IHW'}
            A{n}=logical(A{n}-A{n}.*W{n}); % Remove missing links
        case {'IWRM','IWDM','IWHW'}
            A{n}=A{n}-A{n}.*W{n}; % Remove missing links
    end
    [Iw{n},Jw{n}]=find(W{n});
    West{n}=sparse(I,J);
    sumA=sumA+sum(sum(A{n}));
    nnzA=nnzA+nnz(A{n});
end

% Initialize Z
ind=ceil(noc*rand(1,I));
ind(1:noc) = 1:noc;
Z=mgetopt(opts,'Z',accumarray([ind' (1:J)'],ones(J,1),[noc J]));
noc=size(Z,1);

% Zet remaining parameters
alpha=mgetopt(opts,'alpha',log(J));
verbose=mgetopt(opts,'verbose',1);
switch method
    case {'IRM','IDM','IHW'}
        eta0p=mgetopt(opts,'eta0p',[1 1]); % pseudo counts of links between clusters
        eta0n=mgetopt(opts,'eta0n',[1 1]); % pseudo counts of non-links between clusters
    case {'IWRM','IWDM','IWHW'}
        eta0p=mgetopt(opts,'eta0p',[1 1]); % pseudo counts of links between clusters
        eta0n=mgetopt(opts,'eta0n',[1 1]); % pseudo counts of total entries between clusters
end
dZstep=mgetopt(opts,'dZstep',25); % pseudo counts of links between clusters
init_sample_iter=mgetopt(opts,'init_sample_iter',10);
nsampleiter=mgetopt(opts,'nsampleiter',10);
maxiter=init_sample_iter+nsampleiter;

sample=struct;
L=zeros(1,maxiter);
cpu_time=L;
sstep=0;
westiter=0;
vv1=cell(1,N);
for n=1:N
    [ii,jj,vv1{n}]=find(W{n}.*A{n}+W{n});    
    vv1{n}=vv1{n}-1;
end

iter=0;
if verbose % Display algorithm    
    disp(['Non-parametric clustering based on the Infinite ' method ' for ' type ' graphs'])
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
        
    % Gibbs sampling of Z in random order    
    [Z,logP_A,logP_Z]=Gibbs_sample_ZIRM(Z,A,W,eta0p,eta0n,alpha,randperm(J),type,method);            
    % Zplit/merge sample Z
    %for rep=1:noc
         [Z,logP_A,logP_Z]=split_merge_sample_Z(Z,A,W,eta0p,eta0n,alpha,logP_A,logP_Z,type,method);
    %end     
    noc=size(Z,1);        
    
    eta=estimateRho(A,W,Z,eta0p,eta0n,method,type);
    
    Q=logP_A+logP_Z;
    dQ=Q-Qold;    
    L(iter)=Q;
    t_iter=toc;
    cpu_time(iter)=t_iter;    
    if mod(iter,dZstep)==0 && iter>=init_sample_iter
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
    if  iter>init_sample_iter 
        westiter=westiter+1;        
        if iter==maxiter-nsampleiter+1
            disp(['Initiating estimation of missing links for the last ' num2str(nsampleiter) ' iteration(s)']);   
        end
        step=10000;
        for n=1:N
            val=zeros(1,length(Iw{n}));
            for k=1:ceil((length(Iw{n})/step))
                 ind=(k-1)*step+1:min([k*step, length(Iw{n})]);   
                val(ind)=sum(Z(:,Iw{n}(ind)).*(eta(:,:,n)*Z(:,Jw{n}(ind))))+eps;
            end
            West{n}=West{n}+sparse(Iw{n},Jw{n},val,I,J)/nsampleiter;                
            if method(2)=='W'
                predL(n)=predL(n)+sum(vv1{n}.*log(val')-val')/(nsampleiter*length(val));
            else                
                predL(n)=predL(n)+sum(vv1{n}.*log(val')+(1-vv1{n}).*log(1-val'))/(nsampleiter*length(val));
            end
        end
    end
        
end
% sort Z
[val,ind]=sort(sum(Z,2),'descend');
Z=Z(ind,:);
eta=eta(ind,:);
eta=eta(:,ind);
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
function eta=estimateRho(A,W,Z,eta0p,eta0n,method,type)
        
        sumZ=sum(Z,2);
        [noc, I]=size(Z);
        N=length(A);
        Ap=eta0p(1)*eye(noc)+eta0p(2)*(ones(noc)-eye(noc));
        An=eta0n(1)*eye(noc)+eta0n(2)*(ones(noc)-eye(noc));    
        eta=zeros(noc,noc,N);
        ZZT=sumZ*sumZ'-diag(sumZ);
        ZZT2=sumZ.*(sumZ-1);
        for n=1:N            
                switch method
                   case 'IRM'
                        ZAZt=Z*A{n}*Z';
                        ZWZt=Z*W{n}*Z';
                        if strcmp(type,'UnDirected')                    
                            n_link=ZAZt+ZAZt';    
                            n_link=n_link-0.5*diag(diag(n_link));
                            n_link=n_link+Ap;
                            n_nonlink=ZZT-ZAZt-ZWZt-ZAZt'-ZWZt';
                            n_nonlink=n_nonlink-0.5*diag(diag(n_nonlink));
                            n_nonlink=n_nonlink+An;                   
                        else                    
                            n_link=ZAZt'+Ap;                                            
                            n_nonlink=ZZT-ZAZt-ZWZt+An;                                        
                        end
                        eta(:,:,n)=n_link./(n_link+n_nonlink);
                    case 'IWRM'
                        ZAZt=Z*A{n}*Z';
                        ZWZt=Z*W{n}*Z';
                        if strcmp(type,'UnDirected')                    
                            n_link=ZAZt+ZAZt';    
                            n_link=n_link-0.5*diag(diag(n_link));
                            n_link=n_link+Ap;
                            n_nonlink=ZZT-ZWZt-ZWZt';
                            n_nonlink=n_nonlink-0.5*diag(diag(n_nonlink));
                            n_nonlink=n_nonlink+An;                   
                        else                    
                            n_link=ZAZt'+Ap;                                            
                            n_nonlink=ZZT-ZWZt+An;                                        
                        end
                        eta(:,:,n)=n_link./(n_nonlink);
                    case 'IDM'
                        ZAZt=sum((Z*A{n}).*Z,2);
                        ZWZt=sum((Z*W{n}).*Z,2);
                        if strcmp(type,'UnDirected')                            
                            eta(:,:,n)=diag((ZAZt+eta0p(1))./(ZZT2/2-ZWZt+eta0p(1)+eta0n(1)))+(sum(sum(A{n}))-sum(ZAZt)+eta0p(2))/(I*(I-1)/2-sum(ZZT2/2)-(nnz(W{n})-sum(ZWZt))+eta0p(2)+eta0n(2))*(ones(noc)-eye(noc));
                        else
                            eta(:,:,n)=diag((ZAZt+eta0p(1))./(ZZT2-ZWZt+eta0p(1)+eta0n(1)))+(sum(sum(A{n}))-sum(ZAZt)+eta0p(2))/(I*(I-1)-sum(ZZT2)-(nnz(W{n})-sum(ZWZt))+eta0p(2)+eta0n(2))*(ones(noc)-eye(noc));
                        end
                     case 'IWDM'
                        ZAZt=sum((Z*A{n}).*Z,2);                        
                        ZWZt=sum((Z*W{n}).*Z,2);
                        if strcmp(type,'UnDirected')                            
                            eta(:,:,n)=diag((ZAZt+eta0p(1))./(ZZT2/2-ZWZt+eta0n(1)))+(sum(sum(A{n}))-sum(ZAZt)+eta0p(2))/(I*(I-1)/2-sum(ZZT2/2)-(nnz(W{n})-sum(ZWZt))+eta0n(2))*(ones(noc)-eye(noc));
                        else
                            eta(:,:,n)=diag((ZAZt+eta0p(1))./(ZZT2-ZWZt+eta0n(1)))+(sum(sum(A{n}))-sum(ZAZt)+eta0p(2))/(I*(I-1)-sum(ZZT2)-(nnz(W{n})-sum(ZWZt))+eta0n(2))*(ones(noc)-eye(noc));
                        end
                    case 'IHW'
                        ZAZt=sum((Z*A{n}).*Z,2);
                        ZWZt=sum((Z*W{n}).*Z,2);
                        if strcmp(type,'UnDirected')                            
                            eta(:,:,n)=(sum(ZAZt)+eta0p(1))./(sum(ZZT2/2-ZWZt)+eta0p(1)+eta0n(1))*eye(noc)+(sum(sum(A{n}))-sum(ZAZt)+eta0p(2))/(I*(I-1)/2-sum(ZZT2/2)-(nnz(W{n})-sum(ZWZt))+eta0p(2)+eta0n(2))*(ones(noc)-eye(noc));
                        else
                            eta(:,:,n)=(sum(ZAZt)+eta0p(1))./(sum(ZZT2-ZWZt)+eta0p(1)+eta0n(1))*eye(noc)+(sum(sum(A{n}))-sum(ZAZt)+eta0p(2))/(I*(I-1)-sum(ZZT2)-(nnz(W{n})-sum(ZWZt))+eta0p(2)+eta0n(2))*(ones(noc)-eye(noc));
                        end
                    case 'IWHW'
                        ZAZt=sum((Z*A{n}).*Z,2);                        
                        ZWZt=sum((Z*W{n}).*Z,2);
                        if strcmp(type,'UnDirected')                            
                            eta(:,:,n)=(sum(ZAZt)+eta0p(1))/(sum(ZZT2/2-ZWZt)+eta0n(1))*eye(noc)+(sum(sum(A{n}))-sum(ZAZt)+eta0p(2))/(I*(I-1)/2-sum(ZZT2/2)-(nnz(W{n})-sum(ZWZt))+eta0n(2))*(ones(noc)-eye(noc));
                        else
                            eta(:,:,n)=(sum(ZAZt)+eta0p(1))/(sum(ZZT2-ZWZt)+eta0n(1))*eye(noc)+(sum(sum(A{n}))-sum(ZAZt)+eta0p(2))/(I*(I-1)-sum(ZZT2)-(nnz(W{n})-sum(ZWZt))+eta0n(2))*(ones(noc)-eye(noc));
                        end
                end           
        end
        
% -------------------------------------------------------------------------  
function [Z,logP_A,logP_Z]=split_merge_sample_Z(Z,A,W,eta0p,eta0n,alpha,logP_A,logP_Z,type,method)

    
    %[logP_A_t,logP_Z_t]=evalLikelihood(Z,A,W,eta0p,eta0n,alpha,type,method);   
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
                [Z_t,logP_A_t,logP_Z_t,logQ_trans,comp]=Gibbs_sample_ZIRM(Z_t,A,W,eta0p,eta0n,alpha,setZ(randperm(n_setZ)),type,method,comp);                        
            end     
        else
           logQ_trans=0;
           [logP_A_t,logP_Z_t]=evalLikelihood(Z_t,A,W,eta0p,eta0n,alpha,type,method);                 
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
        [logP_A_t,logP_Z_t]=evalLikelihood(Z_t,A,W,eta0p,eta0n,alpha,type,method);                

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
                [Z_tt,logP_A_tt,logP_Z_tt,logQ_trans,comp]=Gibbs_sample_ZIRM(Z_tt,A,W,eta0p,eta0n,alpha,setZ(randperm(n_setZ)),type,method,comp);               
            end
            Force=[1 2]*Z([clust1 clust2],:);        
            [Z_tt,logP_A_tt,logP_Z_tt,logQ_trans]=Gibbs_sample_ZIRM(Z_tt,A,W,eta0p,eta0n,alpha,setZ(randperm(n_setZ)),type,method,comp,Force);                        
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


% -------------------------------------------------------------------------  
function [logP_A,logP_Z]=evalLikelihood(Z,A,W,eta0p,eta0n,alpha,type,method)

    N=length(A);
    [I,J]=size(A{1});
    noc=size(Z,1);
    sumZ=sum(Z,2);        
    logP_A=0;
    ii=find(triu(ones(noc)));      
  
    ZZT=sumZ*sumZ'-diag(sumZ);
    ZZT2=sumZ.*(sumZ-1);
    for n=1:N        
        if strcmp(type,'UnDirected')            
            switch method
                case 'IRM'
                    ZAZt=Z*A{n}*Z';
                    ZWZt=Z*W{n}*Z';      
                    n_link=ZAZt+ZAZt';    
                    n_link=n_link-0.5*diag(diag(n_link));
                    Ap=eta0p(1)*eye(noc)+eta0p(2)*(ones(noc)-eye(noc));
                    An=eta0n(1)*eye(noc)+eta0n(2)*(ones(noc)-eye(noc));       
                    n_link=n_link+Ap;
                    n_nonlink=ZZT-ZAZt-ZWZt-ZAZt'-ZWZt';
                    n_nonlink=n_nonlink-0.5*diag(diag(n_nonlink));
                    n_nonlink=n_nonlink+An;
                    logP_A=logP_A+sum(betaln(n_link(ii),n_nonlink(ii)))-sum([noc noc*(noc-1)/2].*betaln(eta0p,eta0n));
                case 'IDM'
                    ZAZt=sum((Z*A{n}).*Z,2);
                    ZWZt=sum((Z*W{n}).*Z,2);                          
                    n_link=[ZAZt+eta0p(1); sum(sum(A{n}))-sum(ZAZt)+eta0p(2)];                    
                    n_nonlink=[ZZT2/2-ZAZt-ZWZt+eta0n(1); I*(I-1)/2-sum(ZZT2)/2-(sum(sum(A{n}))-sum(ZAZt))-(nnz(W{n})-(sum(ZWZt)))+eta0n(2)];                    
                    logP_A=logP_A+sum(betaln(n_link,n_nonlink))-sum([noc 1].*betaln(eta0p,eta0n));
                case 'IHW'
                    ZZT2=sum(ZZT2);
                    ZAZt=sum(sum((Z*A{n}).*Z,2));
                    ZWZt=sum(sum((Z*W{n}).*Z,2));                          
                    n_link=[ZAZt+eta0p(1); sum(sum(A{n}))-ZAZt+eta0p(2)];                    
                    n_nonlink=[ZZT2/2-ZAZt-ZWZt+eta0n(1); I*(I-1)/2-ZZT2/2-(sum(sum(A{n}))-ZAZt)-(nnz(W{n})-ZWZt)+eta0n(2)];                    
                    logP_A=logP_A+sum(betaln(n_link,n_nonlink))-sum(betaln(eta0p,eta0n));
                case 'IWRM'
                    ZAZt=Z*A{n}*Z';
                    ZWZt=Z*W{n}*Z';      
                    n_link=ZAZt+ZAZt';    
                    n_link=n_link-0.5*diag(diag(n_link));
                    Ap=eta0p(1)*eye(noc)+eta0p(2)*(ones(noc)-eye(noc));
                    An=eta0n(1)*eye(noc)+eta0n(2)*(ones(noc)-eye(noc));       
                    n_link=n_link+Ap;
                    n_nonlink=ZZT-ZWZt-ZWZt';
                    n_nonlink=n_nonlink-0.5*diag(diag(n_nonlink));
                    n_nonlink=n_nonlink+An;
                    logP_A=logP_A+sum(clusterEvalPoisson(n_link(ii),n_nonlink(ii)))-sum([noc noc*(noc-1)/2].*clusterEvalPoisson(eta0p,eta0n));
                case 'IWDM'
                    ZAZt=sum((Z*A{n}).*Z,2);
                    ZWZt=sum((Z*W{n}).*Z,2);                          
                    n_link=[ZAZt+eta0p(1); sum(sum(A{n}))-sum(ZAZt)+eta0p(2)];                    
                    n_nonlink=[ZZT2/2-ZWZt+eta0n(1); I*(I-1)/2-sum(ZZT2)/2-(nnz(W{n})-(sum(ZWZt)))+eta0n(2)];                    
                    logP_A=logP_A+sum(clusterEvalPoisson(n_link,n_nonlink))-sum([noc 1].*clusterEvalPoisson(eta0p,eta0n));
                case 'IWHW'
                    ZZT2=sum(ZZT2);
                    ZAZt=sum(sum((Z*A{n}).*Z,2));
                    ZWZt=sum(sum((Z*W{n}).*Z,2));                          
                    n_link=[ZAZt+eta0p(1); sum(sum(A{n}))-ZAZt+eta0p(2)];                    
                    n_nonlink=[ZZT2/2-ZWZt+eta0n(1); I*(I-1)/2-ZZT2/2-(nnz(W{n})-ZWZt)+eta0n(2)];                    
                    logP_A=logP_A+sum(clusterEvalPoisson(n_link,n_nonlink))-sum(clusterEvalPoisson(eta0p,eta0n));
            end            
        else            
            switch method
                case 'IRM'
                    ZAZt=Z*A{n}*Z';
                    ZWZt=Z*W{n}*Z';
                    Ap=eta0p(1)*eye(noc)+eta0p(2)*(ones(noc)-eye(noc));
                    An=eta0n(1)*eye(noc)+eta0n(2)*(ones(noc)-eye(noc));       
                    n_link=ZAZt+Ap;
                    n_nonlink=ZZT-ZAZt-ZWZt+An;                                                             
                    logP_A=logP_A+sum(sum(betaln(n_link,n_nonlink)))-sum([noc noc*(noc-1)].*betaln(eta0p,eta0n));
                case 'IDM'
                    ZAZt=sum((Z*A{n}).*Z,2);
                    ZWZt=sum((Z*W{n}).*Z,2);                          
                    n_link=[ZAZt+eta0p(1); sum(sum(A{n}))-sum(ZAZt)+eta0p(2)];                    
                    n_nonlink=[ZZT2-ZAZt-ZWZt+eta0n(1); I*(I-1)-sum(ZZT2)-(sum(sum(A{n}))-sum(ZAZt))-(nnz(W{n})-(sum(ZWZt)))+eta0n(2)];                    
                    logP_A=logP_A+sum(betaln(n_link,n_nonlink))-sum([noc 1].*betaln(eta0p,eta0n));
                case 'IHW'
                    ZZT2=sum(ZZT2);
                    ZAZt=sum(sum((Z*A{n}).*Z,2));
                    ZWZt=sum(sum((Z*W{n}).*Z,2));                          
                    n_link=[ZAZt+eta0p(1); sum(sum(A{n}))-ZAZt+eta0p(2)];                    
                    n_nonlink=[ZZT2-ZAZt-ZWZt+eta0n(1); I*(I-1)-ZZT2-(sum(sum(A{n}))-ZAZt)-(nnz(W{n})-ZWZt)+eta0n(2)];                    
                    logP_A=logP_A+sum(betaln(n_link,n_nonlink))-sum(betaln(eta0p,eta0n));
                case 'IWRM'
                    ZAZt=Z*A{n}*Z';
                    ZWZt=Z*A{n}*Z';
                    Ap=eta0p(1)*eye(noc)+eta0p(2)*(ones(noc)-eye(noc));
                    An=eta0n(1)*eye(noc)+eta0n(2)*(ones(noc)-eye(noc));       
                    n_link=ZAZt+Ap;
                    n_nonlink=ZZT-ZWZt+An;                                      
                    logP_A=logP_A+sum(sum(clusterEvalPoisson(n_link,n_nonlink)))-sum([noc noc*(noc-1)].*clusterEvalPoisson(eta0p,eta0n));
                case 'IWDM'
                    ZAZt=sum((Z*A{n}).*Z,2);
                    ZWZt=sum((Z*W{n}).*Z,2);                          
                    n_link=[ZAZt+eta0p(1); sum(sum(A{n}))-sum(ZAZt)+eta0p(2)];                    
                    n_nonlink=[ZZT2-ZWZt+eta0n(1); I*(I-1)-sum(ZZT2)-(nnz(W{n})-(sum(ZWZt)))+eta0n(2)];                    
                    logP_A=logP_A+sum(clusterEvalPoisson(n_link,n_nonlink))-sum([noc 1].*clusterEvalPoisson(eta0p,eta0n));          
                case 'IWHW'                    
                    ZZT2=sum(ZZT2);
                    ZAZt=sum(sum((Z*A{n}).*Z,2));
                    ZWZt=sum(sum((Z*W{n}).*Z,2));                          
                    n_link=[ZAZt+eta0p(1); sum(sum(A{n}))-ZAZt+eta0p(2)];                    
                    n_nonlink=[ZZT2-ZWZt+eta0n(1); I*(I-1)-ZZT2-(nnz(W{n})-ZWZt)+eta0n(2)];                    
                    logP_A=logP_A+sum(clusterEvalPoisson(n_link,n_nonlink))-sum(clusterEvalPoisson(eta0p,eta0n));          
            end
        end
        if sum(sum(betaln(n_link,n_nonlink)))==Inf
            keyboard
        end
    end
    logP_Z=noc*log(alpha)+sum(gammaln(full(sumZ)))-gammaln(J+alpha)+gammaln(alpha);    
    
% -------------------------------------------------------------------------    
function  B=clusterEvalPoisson(Np,Nn)

    B=gammaln(Np)-Np.*log(Nn);

% -------------------------------------------------------------------------
function [Z,logP_A,logP_Z,logQ_trans,comp]=Gibbs_sample_ZIRM(Z,A,W,eta0p,eta0n,alpha,JJ,type,method,comp,Force)        
    if nargin<12
        Force=[];
    end
    if nargin<11
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
            
    Ap=eta0p(1)*eye(noc)+eta0p(2)*(ones(noc)-eye(noc));
    An=eta0n(1)*eye(noc)+eta0n(2)*(ones(noc)-eye(noc));    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Count initial number of links and non-links between groups
    switch method
        case {'IRM','IWRM'}
            n_link=zeros(noc,noc,N);
            n_nonlink=n_link;
        case {'IDM','IWDM','IWHW','IHW'}
            n_link=zeros(noc+1,N);
            n_nonlink=n_link;            
    end
    if ~strcmp(type,'UnDirected')
        A2=cell(1,N);
        W2=cell(1,N);
    end
    ZZT=sumZ*sumZ'-diag(sumZ);
    ZZT2=sumZ.*(sumZ-1);
    for n=1:N        
        if strcmp(type,'UnDirected')
            switch method
                case 'IRM'   
                    ZAZt=Z*A{n}*Z';
                    ZWZt=Z*W{n}*Z';
                    n_link(:,:,n)=ZAZt+ZAZt';    
                    n_link(:,:,n)=n_link(:,:,n)-0.5*diag(diag(n_link(:,:,n)));
                    n_link(:,:,n)=n_link(:,:,n)+Ap;
                    n_nonlink(:,:,n)=ZZT-ZAZt-ZWZt-ZAZt'-ZWZt';
                    n_nonlink(:,:,n)=n_nonlink(:,:,n)-0.5*diag(diag(n_nonlink(:,:,n)));
                    n_nonlink(:,:,n)=n_nonlink(:,:,n)+An;        
                case 'IDM'
                    ZAZt=sum((Z*A{n}).*Z,2);
                    ZWZt=sum((Z*W{n}).*Z,2);
                    n_link(:,n)=[ZAZt+eta0p(1); (nnz(A{1})-sum(ZAZt))+eta0p(2)];                                                                
                    nnzC=ZZT2/2;
                    T=nnzC-ZAZt-ZWZt;                                        
                    n_nonlink(:,n)=[T+eta0n(1); I*(I-1)/2-sum(nnzC)-(sum(sum(A{n}))-sum(ZAZt))-(nnz(W{n})-sum(ZWZt))+eta0n(2)];                                                                                                        
                case 'IHW'
                    ZAZt=sum((Z*A{n}).*Z,2);
                    ZWZt=sum((Z*W{n}).*Z,2);
                    n_link(:,n)=[ZAZt; (nnz(A{1})-sum(ZAZt))];                                                                
                    nnzC=ZZT2/2;
                    T=nnzC-ZAZt-ZWZt;                                        
                    n_nonlink(:,n)=[T; I*(I-1)/2-sum(nnzC)-(sum(sum(A{n}))-sum(ZAZt))-(nnz(W{n})-sum(ZWZt))];                                                                                                        
                case 'IWRM'
                    ZAZt=Z*A{n}*Z';
                    ZWZt=Z*W{n}*Z';
                    n_link(:,:,n)=ZAZt+ZAZt';    
                    n_link(:,:,n)=n_link(:,:,n)-0.5*diag(diag(n_link(:,:,n)));
                    n_link(:,:,n)=n_link(:,:,n)+Ap;
                    n_nonlink(:,:,n)=ZZT-ZWZt-ZWZt';
                    n_nonlink(:,:,n)=n_nonlink(:,:,n)-0.5*diag(diag(n_nonlink(:,:,n)));
                    n_nonlink(:,:,n)=n_nonlink(:,:,n)+An;        
                case 'IWDM'
                    ZAZt=sum((Z*A{n}).*Z,2);
                    ZWZt=sum((Z*W{n}).*Z,2);
                    n_link(:,n)=[ZAZt+eta0p(1); sum(sum(A{n}))-sum(ZAZt)+eta0p(2)];                                                                                                                               
                    n_nonlink(:,n)=[ZZT2/2-ZWZt+eta0n(1); I*(I-1)/2-sum(ZZT2)/2-(nnz(W{n})-sum(ZWZt))+eta0n(2)];                                                                                                                           
                case 'IWHW'
                    ZAZt=sum((Z*A{n}).*Z,2);
                    ZWZt=sum((Z*W{n}).*Z,2);
                    n_link(:,n)=[ZAZt; sum(sum(A{n}))-sum(ZAZt)];                                                                                                                               
                    n_nonlink(:,n)=[ZZT2/2-ZWZt; I*(I-1)/2-sum(ZZT2)/2-(nnz(W{n})-sum(ZWZt))]; 
            end  
            switch method
                case {'IRM','IDM','IHW'}
                    A{n}=logical(A{n}+A{n}');
                case {'IWRM','IWDM','IWHW'}
                    A{n}=A{n}+A{n}';
            end
            W{n}=logical(W{n}+W{n}');
        else            
            A2{n}=A{n}';
            W2{n}=W{n}';            
            switch method
                case 'IRM'                             
                    ZAZt=Z*A{n}*Z';
                    ZWZt=Z*W{n}*Z';                    
                    n_link(:,:,n)=ZAZt+Ap;    
                    n_nonlink(:,:,n)=ZZT-ZAZt-ZWZt+An;                        
                case 'IDM'
                    ZAZt=sum((Z*A{n}).*Z,2);
                    ZWZt=sum((Z*W{n}).*Z,2);
                    n_link(:,n)=[ZAZt+eta0p(1); (sum(sum(A{n}))-sum(ZAZt))+eta0p(2)];                                                                
                    nnzC=ZZT2;
                    T=nnzC-ZAZt-ZWZt;                                        
                    n_nonlink(:,n)=[T+eta0n(1); I*(I-1)-sum(nnzC)-(sum(sum(A{n}))-sum(ZAZt))-(nnz(W{n})-sum(ZWZt))+eta0n(2)];                                                                                                                           
                case 'IHW'
                    ZAZt=sum((Z*A{n}).*Z,2);
                    ZWZt=sum((Z*W{n}).*Z,2);
                    n_link(:,n)=[ZAZt; (sum(sum(A{n}))-sum(ZAZt))];                                                                
                    nnzC=ZZT2;
                    T=nnzC-ZAZt-ZWZt;                                        
                    n_nonlink(:,n)=[T; I*(I-1)-sum(nnzC)-(sum(sum(A{n}))-sum(ZAZt))-(nnz(W{n})-sum(ZWZt))];                                                                                                                           
                case 'IWRM'
                    ZAZt=Z*A{n}*Z';
                    ZWZt=Z*W{n}*Z';                    
                    n_link(:,:,n)=ZAZt+Ap;    
                    n_nonlink(:,:,n)=ZZT-ZWZt+An;                        
                case 'IWDM'
                    ZAZt=sum((Z*A{n}).*Z,2);
                    ZWZt=sum((Z*W{n}).*Z,2);
                    n_link(:,n)=[ZAZt+eta0p(1); (sum(sum(A{n}))-sum(ZAZt))+eta0p(2)];                                                                                    
                    n_nonlink(:,n)=[ZZT2-ZWZt+eta0n(1); I*(I-1)-sum(ZZT2)-(nnz(W{n})-sum(ZWZt))+eta0n(2)];                                                                                                        
                case 'IWHW'
                    ZAZt=sum((Z*A{n}).*Z,2);
                    ZWZt=sum((Z*W{n}).*Z,2);
                    n_link(:,n)=[ZAZt; (sum(sum(A{n}))-sum(ZAZt))];                                                                                    
                    n_nonlink(:,n)=[ZZT2-ZWZt; I*(I-1)-sum(ZZT2)-(nnz(W{n})-sum(ZWZt))];                                                                                                        
            end                
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

                        e_noc=ones(1,noc);
                        link1=zeros(noc,noc+1,N);
                        nolink1=zeros(noc,noc+1,N);
                        link2=zeros(noc,noc+1,N);
                        nolink2=zeros(noc,noc+1,N);                
                        for n=1:N
                            link1(:,1:noc,n)=n_link(:,:,n)+ZA1k(:,n*e_noc)+diag(ZA2k(:,n));                             
                            nolink1(:,1:noc,n)=n_nonlink(:,:,n)+nZA1k(:,n*e_noc)+diag(nZA2k(:,n));
                            link2(:,1:noc,n)=n_link(:,:,n)'+ZA2k(:,n*e_noc)+diag(ZA1k(:,n));
                            nolink2(:,1:noc,n)=n_nonlink(:,:,n)'+nZA2k(:,n*e_noc)+diag(nZA1k(:,n));
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
                            link(1:noc,1,:)=link_within(enoc',:)+ZA1k+eta0p(1);
                            link(1:noc,2,:)=n_link((noc+1)*enoc',:)+enoc'*sum(ZA1k,1)-ZA1k+eta0p(2);
                            link(noc+1,1,:)=link_within+eta0p(1);                                                        
                            link(noc+1,2,:)=n_link(noc+1,:)+sum(ZA1k,1)+eta0p(2);
                            nolink(1:noc,1,:)=nonlink_within(enoc',:)+nZA1k+eta0n(1);
                            nolink(1:noc,2,:)=n_nonlink((noc+1)*enoc',:)+enoc'*sum(nZA1k,1)-nZA1k+eta0n(2);
                            nolink(noc+1,1,:)=nonlink_within+eta0p(1);                            
                            nolink(noc+1,2,:)=n_nonlink(noc+1,:)+sum(nZA1k,1)+eta0p(2);
                        else
                            link(1:noc,1,:)=link_within(enoc',:)+ZA1k+ZA2k+eta0p(1);
                            link(1:noc,2,:)=n_link((noc+1)*enoc',:)+enoc'*sum(ZA1k,1)-ZA1k+enoc'*sum(ZA2k,1)-ZA2k+eta0p(2);
                            link(noc+1,1,:)=link_within+eta0p(1);
                            link(noc+1,2,:)=n_link((noc+1),:)+sum(ZA1k,1)+sum(ZA2k,1)+eta0p(2);
                            nolink(1:noc,1,:)=nonlink_within(enoc',:)+nZA1k+nZA2k+eta0n(1);
                            nolink(1:noc,2,:)=n_nonlink((noc+1)*enoc',:)+enoc'*sum(nZA1k,1)-nZA1k+enoc'*sum(nZA2k,1)-nZA2k+eta0n(2);
                            nolink(noc+1,1,:)=nonlink_within+eta0n(1);
                            nolink(noc+1,2,:)=n_nonlink((noc+1),:)+sum(nZA1k,1)+sum(nZA2k,1)+eta0n(2);
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
                        cluster_eval_d1(noc,:)=0;                                
                        ZA1k(noc,:)=0;
                        nZA1k(noc,:)=0;              
                        if ~strcmp(type,'UnDirected')
                            if N==1
                                cluster_eval_d2=cluster_eval_d_2(:,noc);
                            else
                                cluster_eval_d2=permute(cluster_eval_d_2(:,noc,:),[1 3 2]);                                                
                            end                    
                            cluster_eval_d2(noc,:)=0;
                            ZA2k(noc,:)=0;
                            nZA2k(noc,:)=0;              
                            logQf=logQ(noc)+N*2*(noc-1)*off_const;
                        else                
                            logQf=logQ(noc)+N*(noc-1)*off_const;
                        end
                    case {'IDM','IWDM'}
                        n_link(noc+1,:)=n_link(noc,:);
                        n_link(noc,:)=eta0p(1);
                        n_nonlink(noc+1,:)=n_nonlink(noc,:);
                        n_nonlink(noc,:)=eta0n(1);                                                
                        cluster_eval(noc,:)=0;                                        
                        ZA1k(noc,:)=0;                            
                        nZA1k(noc,:)=0;              
                        if ~strcmp(type,'UnDirected')                            
                            ZA2k(noc,:)=0;                            
                            nZA2k(noc,:)=0;                  
                        end
                        logQf=logQ(noc);                       
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
                                nolink(t,:,2)=n_nonlink(noc+1,:)+sum(nZA1k,1)-nZA1k(comp(t),:)+sum(nZA1k,1)-nZA1k(comp(t),:);
                            end
                        end
                        logQ=sum(sum(clustFun(link,nolink),3),2); 
                    case {'IHW','IWHW'}
                        link=zeros(2,N,2);
                        nolink=link;
                        e=ones(2,1);
                        link_within=sum(n_link(1:noc,:),1);
                        nonlink_within=sum(n_nonlink(1:noc,:),1);
                        if strcmp(type,'UnDirected')
                            for t=1:2                                
                                link(t,:,1)=link_within(e,:)+ZA1k(comp(t),:)+eta0p(1);
                                link(t,:,2)=n_link(noc+1,:)+sum(ZA1k,1)-ZA1k(comp(t),:)+eta0p(2);
                                nolink(t,:,1)=nonlink_within(e,:)+nZA1k(comp(t),:)+etanp(1);
                                nolink(t,:,2)=n_nonlink(noc+1,:)+sum(nZA1k,1)-nZA1k(comp(t),:)+eta0n(2);
                            end
                        else                                                        
                            for t=1:2
                                link(t,:,1)=link_within(e,:)+ZA1k(comp(t),:)+ZA2k(comp(t),:)+eta0p(1);
                                link(t,:,2)=n_link(noc+1,:)+sum(ZA1k,1)-ZA1k(comp(t),:)+sum(ZA2k,1)-ZA2k(comp(t),:)+eta0p(2);
                                nolink(t,:,1)=nonlink_within(e,:)+nZA1k(comp(t),:)+nZA2k(comp(t),:)+eta0n(1);
                                nolink(t,:,2)=n_nonlink(noc+1,:)+sum(nZA1k,1)-nZA1k(comp(t),:)+sum(nZA1k,1)-nZA1k(comp(t),:)+eta0n(2);
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
            logQ_trans=logQ_trans+log(QQ(ind)+eps);
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
                    n_nonlink(end,:)=n_nonlink(end,:)+sum(nZA1k)-nZA1k(ind,:);                            
                else
                    n_link(ind,:)=n_link(ind,:)+ZA1k(ind,:)+ZA2k(ind,:);                    
                    n_nonlink(ind,:)=n_nonlink(ind,:)+nZA1k(ind,:)+nZA2k(ind,:);                          
                    n_link(end,:)=n_link(end,:)+sum(ZA1k,1)-ZA1k(ind,:)+sum(ZA2k,1)-ZA2k(ind,:);                    
                    n_nonlink(end,:)=n_nonlink(end,:)+sum(nZA1k)-nZA1k(ind,:)+sum(nZA2k,1)-nZA2k(ind,:);                            
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
            logP_A=sum(sum(cluster_eval))+sum(clustFun(n_link(noc+1,:),n_nonlink(noc+1,:)))-N*sum([noc 1].*[diag_const off_const]);                           
        case {'IHW','IWHW'}
            logP_A=logQf-N*sum([diag_const off_const]);                           
    end
        
    