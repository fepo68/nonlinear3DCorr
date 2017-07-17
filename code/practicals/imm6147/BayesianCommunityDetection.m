function [L,cpu_time,Z,eta,gap,sample,West,predL]=BayesianCommunityDetection(A,W,noc,opts)

% Non-parametric clustering of un-directed graphs based on collapsed Gibbs sampling with split
% merge
%
% Usage: 
%  [L,cpu_time,Z,eta,sample,West]=BayesianCommunityDetection(A,W,noc,opts)
%
% Input:
%   A       cell array of I x I  sparse adjacency matrix (can be triuangular upper matrix)
%   W       cell array of I x I  sparse missing value indicator matrix (can be triuangular upper matrix)
%   noc     number of clusters
%   opts.
%           init_sample_iter initial number of sampling iterations  (i.e. burnin)
%           nsampleiter iterations to use for prediction
%           Z           I x noc matrix of community identities
%           dSstep      number of iteratiosn between each recorded sample
%           verbose     1: display iteration results, 0: no display
%           rho0p       1x2 vector of pseudo link counts wihtin and 
%                       between communities (default: [1 1])
%           rho0n       1x2 vector of pseudo non-link counts wihtin and 
%                       between communities and (default: [1 1])
%           gap_type    'same' or 'individual', (default: 'same')
%           constGap    make gap constant (default: false)
%           sameEta     use same value of eta for multiple graphs (default: false)

% Output:
%   L           Log-likelihood function at each iteration
%   cpu_time    cpu-time cost for each iteration
%   Z           Estimated clustering assigment matrix
%   eta         Estimated relations
%   gap         Estimated Gap
%   sample      sampled parameters at each dSstep iterationCluster indicator matrix of size d x I where d<=noc before
%   West        Estimated link probabilities for missing links and non-links
%   predL       Predictive log-likelihood based on the use of the expected
%               value of eta rather than expected value of log(eta) and
%               log(1-eta) given as a sample average, i.e.
%               returns for Bernoulli: \sum_{ij\in W} A_ij*log(eta_ij)+(1-A_ij)*log(1-eta_ij)
%               returns for Poisson: \sum_{ij\in W} A_ij*log(eta_ij)-eta_ij
%
% This code is provided as is without any kind of warranty!
%
%
% Code updates: 
%   12 August 2013:  logQ_trans calculated in the log-domain
%                    deleted unnecessary kode snippets
%                    added diag_const in Gibbs sampler and changed line 505 where 'rho0p(1)-rho0n(2)' changed to 'rho0n(1)-rho0n(2)' this was only an issue if priors 
%                    were set differently than the default values
%   16 Oktober 2013: fixed bug in MHsampleGap for Directed networks as well
%                    as potential division by zero estimating rho_noise_t
%   23 Oktober 2013: when density is too low for betaincomplete and
%                    gammaincomplete functions numeric integration implemented for the estimation of rho_noise_t
%   8 December 2013: in line 687 added the line: gap(gap<eps)=eps; to
%                    ensure numeric stability if gap becomes too small
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
type=mgetopt(opts,'type','UnDirected');
method=mgetopt(opts,'method','Binary');
gap_type=mgetopt(opts,'gap_type','same');

par.constGap=mgetopt(opts,'constGap',false);
par.sameEta=mgetopt(opts,'sameEta',false);
Iw=cell(1,par.N);
Jw=cell(1,par.N);
West=cell(1,par.N);
sumA=0;
nnzA=0;
switch method
    case {'Binary'}
        par.rho0p=mgetopt(opts,'rho0p',[1 1]); % pseudo counts of links between clusters
        par.rho0n=mgetopt(opts,'rho0n',[1 1]); % pseudo counts of non-links between clusters        
    case {'Weighted'}
        par.rho0p=mgetopt(opts,'rho0p',[1 1]); % pseudo counts of links between clusters
        par.rho0n=mgetopt(opts,'rho0n',[1 1]); % pseudo counts of non-links between clusters                        
end


for n=1:par.N    
    if strcmp(type,'UnDirected')
        A{n}=triu(A{n});
        W{n}=triu(W{n});
    end    
    W{n}=logical(W{n});
    [ii,jj,vv1{n}]=find(W{n}.*A{n}+W{n});
    vv1{n}=vv1{n}-1;
    switch method
        case {'Binary'}
            A{n}=logical(A{n}-A{n}.*W{n}); % Remove missing links
        case {'Weighted'}
            A{n}=A{n}-A{n}.*W{n}; % Remove missing links
    end
    [Iw{n},Jw{n}]=find(W{n});
    West{n}=sparse(I,J);    
    sumA=sumA+sum(sum(A{n}));
    nnzA=nnzA+nnz(A{n});
    N=par.N;
    if par.sameEta
        if n==1
           At=A{n};
           Wt=W{n};
        else
            At=At+A{n};
            Wt=Wt+W{n};
        end
        if n==par.N
           clear A W
           A{1}=At;
           W{1}=Wt;
        end
        N=1;
    end
end

if strcmp(gap_type,'same')
    gap=mgetopt(opts,'gap',1);
else
    gap=mgetopt(opts,'gap',ones(noc,N));
end

predL=zeros(1,par.N);
rho_diag=zeros(noc,N);
if isfield(opts,'eta')
    for n=1:N
        rho_diag(:,n)=diag(opts.eta(:,:,n));        
    end
else
    for n=1:N
        switch method
            case 'Binary'
                rho_diag(:,n)=betarnd(par.rho0p(1)*ones(noc,1),par.rho0n(1)*ones(noc,1));        
            case 'Weighted'
                rho_diag(:,n)=gamrnd(par.rho0p(1)*ones(noc,1),1/par.rho0n(1)*ones(noc,1));        
        end
        
    end
end

% Initialize Z
ind=ceil(noc*rand(1,I));
ind(1:noc) = 1:noc;
Z=mgetopt(opts,'Z',accumarray([ind' (1:J)'],ones(J,1),[noc J]));
noc=size(Z,1);

% Set remaining parameters
alpha=mgetopt(opts,'alpha',log(J));
verbose=mgetopt(opts,'verbose',1);
dSstep=mgetopt(opts,'dSstep',25); % pseudo counts of links between clusters
init_sample_iter=mgetopt(opts,'init_sample_iter',50);
nsampleiter=mgetopt(opts,'nsampleiter',50);
maxiter=nsampleiter+init_sample_iter;

sample=struct;
L=nan(1,maxiter);
cpu_time=L;
sstep=0;
westiter=0;
    
iter=0;
if verbose % Display algorithm    
    disp(['Non-parametric clustering based on the Bayesian Community Detection using ' method ' method for ' type ' graphs'])
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
    [Z,rho_diag,n_link,n_nonlink,gap]=Gibbs_sample_SBGC(gap,Z,A,W,rho_diag,par,alpha,randperm(J),type,method,gap_type);            

    % Split/merge sample Z
    %for rep=1:noc
         [Z,rho_diag,n_link,n_nonlink,gap]=split_merge_sample_S(gap,Z,A,W,rho_diag,par,alpha,n_link,n_nonlink,type,method,gap_type);
   % end     
    noc=size(Z,1);
 
    % MH sample rho_diag    
    rho_diag=MHsampleRho_diag(gap,n_link,n_nonlink,rho_diag,type,method,par);
    
    % MH sample gap
    [logP_A,logP_S]=evalLikelihood(gap,n_link,n_nonlink,rho_diag,par,alpha,Z,type,method);        
    if ~par.constGap
        [gap,logP_A]=MHsampleGap(gap,gap_type,n_link,n_nonlink,rho_diag,type,method,par.rho0p,par.rho0n,logP_A);        
    end
        
    % form expected value of rho_noise by numeric integration
    rho_noise=zeros(noc,noc,N);
    if strcmp(gap_type,'same')
        rho_cc=gap*minComDens(rho_diag);      
    else
         rho_cc=minComDens(gap.*rho_diag);      
    end
    if strcmp(type,'UnDirected')
        ii=find(triu(ones(noc),1));
    else
        ii=find(ones(noc)-eye(noc));
    end
          
    for n=1:N
        n_link_t=n_link(:,:,n);
        n_nonlink_t=n_nonlink(:,:,n);
        rho_cc_t=rho_cc(:,:,n);
        rho_noise_t=zeros(noc);
        switch method            
            case 'Binary'
                %rho_noise=(beta(n_link+1,n_nonlink).*betainc(rho_cc,n_link+1,n_nonlink))./(beta(n_link,n_nonlink).*betainc(rho_cc,n_link,n_nonlink));                
                rho_noise_t(ii)=n_link_t(ii)./(n_link_t(ii)+n_nonlink_t(ii)).*betainc(rho_cc_t(ii),n_link_t(ii)+1,n_nonlink_t(ii))./(betainc(rho_cc_t(ii),n_link_t(ii),n_nonlink_t(ii))+1e-320);            
                idx=find(rho_noise_t(ii)==0);
                % if density too low use numerical integration
                rho_noise_t(ii(idx))=my_noise_est_beta(rho_cc_t(ii(idx)),n_link_t(ii(idx)),n_nonlink_t(ii(idx)));                                     
            case 'Weighted'
                %rho_noise=(gamma(n_link+1)./(n_nonlink.^(n_link+1)).*gammainc(rho_cc.*n_nonlink,n_link+1))./(gamma(n_link)./(n_nonlink.^n_link).*gammainc(rho_cc.*n_nonlink,n_link));            
                rho_noise_t(ii)=(n_link_t(ii)./n_nonlink_t(ii)).*gammainc(rho_cc_t(ii).*n_nonlink_t(ii),n_link_t(ii)+1)./(gammainc(rho_cc_t(ii).*n_nonlink_t(ii),n_link_t(ii))+1e-320);                            
                idx=find(rho_noise_t(ii)==0);
                % if density is too low use numerical integration
                rho_noise_t(ii(idx))=my_noise_est_gamma(rho_cc_t(ii(idx)),n_link_t(ii(idx)),n_nonlink_t(ii(idx)));                                 
        end
        rho_noise(:,:,n)=rho_noise_t;
    end
           
    eta=zeros(noc,noc,N);
    for n=1:length(A) 
         if strcmp(type,'UnDirected')
             eta(:,:,n)=rho_noise(:,:,n)+rho_noise(:,:,n)'+diag(rho_diag(:,n));         
         else
             eta(:,:,n)=rho_noise(:,:,n)+diag(rho_diag(:,n));         
         end
    end              
    
    Q=logP_A+logP_S;
    dQ=Q-Qold;    
    L(iter)=Q;
    t_iter=toc;
    cpu_time(iter)=t_iter;    
    if mod(iter,dSstep)==0 && iter>=init_sample_iter
        sstep=sstep+1;
        sample.iteration(sstep)=iter;
        sample.Z{sstep}=Z;        
        sample.eta{sstep}=eta;  
        sample.gap{sstep}=gap;  
    end
    if Q>Qbest
        Qbest=Q;
        sample.MAP.L=Q;
        sample.MAP.iteration=iter;
        sample.MAP.Z=Z;
        sample.MAP.eta=eta;
        sample.MAP.gap=gap;
    end
    if rem(iter,1)==0 && verbose        
        disp(sprintf('%12.0f | %12.4e | %12.4e | %12.0f |%12.4f ',iter,Q,dQ/abs(Q),noc,t_iter));
    end
        
    % Estimate missing link probability of sample    
    if  iter>init_sample_iter 
        westiter=westiter+1;        
        if iter==maxiter-nsampleiter+1
            disp(['Initiating estimation of missing links for the last ' num2str(nsampleiter) ' iterations']);   
        end
        step=10000;        
        for n=1:par.N            
            val=zeros(1,length(Iw{n}));
            for k=1:ceil((length(Iw{n})/step))
                ind=(k-1)*step+1:min([k*step, length(Iw{n})]);
                if par.sameEta
                    val(ind)=sum(Z(:,Iw{n}(ind)).*(eta*Z(:,Jw{n}(ind))))+eps;
                else
                    val(ind)=sum(Z(:,Iw{n}(ind)).*(eta(:,:,n)*Z(:,Jw{n}(ind))))+eps;
                end
            end
            West{n}=West{n}+sparse(Iw{n},Jw{n},val,I,J)/nsampleiter;                
            if strcmp(method,'Weighted')
                predL(n)=predL(n)+sum(vv1{n}.*log(val')-val')/(nsampleiter*length(val));
            else                
                predL(n)=predL(n)+sum(vv1{n}.*log(val')+(1-vv1{n}).*log(1-val'))/(nsampleiter*length(val));
            end
        end
    end
     %   figure(1); subplot(1,2,1); plot(L); subplot(1,2,2); imagesc(rho_diag); colorbar;
end

% sort the communities
[val,ind]=sort(sum(Z,2),'descend');
Z=Z(ind,:);
eta=eta(ind,:,:);
eta=eta(:,ind,:);
if ~strcmp(gap_type,'same')
   gap=gap(ind); 
end
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

% -------------------------------------------------------------------------
 function rho_diag=MHsampleRho_diag(gap,n_link,n_nonlink,rho_diag,type,method,par)
        
     rho0p=par.rho0p;
     rho0n=par.rho0n;
     if strcmp(method,'Binary')
        clusterFun=@my_betaub;
        logLike=@lnbetalike;
        drawsample=@betarnd;
     else
        clusterFun=@my_gammaub;
        logLike=@lngammalike;
        drawsample=@gamrnd;
     end
    M=100;
    N=size(rho_diag,2); 
    noc=size(n_link,1);
    accept=0;
    sigma_sq=0.01;
    for rep=1:10  % sample for diagonal elements of eta many times since this is inexpensive
        for n=1:N
            rho_diag_t=rho_diag(:,n);            
            n_link_t=n_link(:,:,n);
            n_nonlink_t=n_nonlink(:,:,n);
            if size(gap,2)>1
               gap_t=gap(:,n); 
            else
               gap_t=gap;
            end
            for d=randperm(noc)
                if strcmp(type,'UnDirected')
                   B=clusterFun(gap_t,n_link_t(:,d),n_nonlink_t(:,d),rho_diag_t,rho0p,rho0n,d);
                else
                   B1=clusterFun(gap_t,n_link_t(:,d),n_nonlink_t(:,d),rho_diag_t,rho0p,rho0n,d);
                   B2=clusterFun(gap_t,n_link_t(d,:)',n_nonlink_t(d,:)',rho_diag_t,rho0p,rho0n,d);
                   B=B1+B2;
                   B(d)=B1(d);
                end
                if rem(rep,3)==1
                    rho_diag_t(d)=drawsample(1,1);
                    deltalogQ=0;
                elseif rem(rep,3)==2
                    a=n_link(d,d);
                    b=n_nonlink(d,d);
                    if strcmp(method,'Weighted')                        
                        rho_diag_t(d)=drawsample(a,1/b);
                    else
                        rho_diag_t(d)=drawsample(a,b);
                    end                    
                    deltalogQ=logLike(rho_diag(d,n),a,b)-logLike(rho_diag_t(d),a,b);                        
                elseif rem(rep,3)==0                                        
                    if strcmp(method,'Weighted')                        
                        mu=rho_diag(d,n);
                        a=mu^2/sigma_sq;
                        b=sigma_sq/mu;            
                        rho_diag_t(d)=drawsample(a,b);
                        mu_new=rho_diag_t(d);
                        a_new=mu_new^2/sigma_sq;
                        b_new=sigma_sq/mu_new;            
                    else                        
                        a=M*rho_diag(d,n);
                        b=M-a;
                        rho_diag_t(d)=drawsample(a,b);
                        a_new=M*rho_diag_t(d);
                        b_new=M-a_new;
                    end                    
                    deltalogQ=logLike(rho_diag(d,n),a_new,b_new)-logLike(rho_diag_t(d),a,b);                        
                end          
                
                if strcmp(type,'UnDirected')
                   Bnew=clusterFun(gap_t,n_link_t(:,d),n_nonlink_t(:,d),rho_diag_t,rho0p,rho0n,d);
                else
                   B1=clusterFun(gap_t,n_link_t(:,d),n_nonlink_t(:,d),rho_diag_t,rho0p,rho0n,d);
                   B2=clusterFun(gap_t,n_link_t(d,:)',n_nonlink_t(d,:)',rho_diag_t,rho0p,rho0n,d);
                   Bnew=B1+B2;
                   Bnew(d)=B1(d); % subtract contribution from d to d which has been included twice (fixed 12 August 2013)
                end
                if rand<exp(sum(Bnew)-sum(B)+deltalogQ) % notice log proposal = 0 hence ignored in MH move
                    accept=accept+1;
                    rho_diag(d,n)=rho_diag_t(d);
                else
                    rho_diag_t(d)=rho_diag(d,n);
                end                
            end
        end
    end
    disp([' Accepted ' num2str(accept) ' samples for eta']);
        

% -------------------------------------------------------------------------  
function [Z,rho_diag,n_link,n_nonlink,gap]=split_merge_sample_S(gap,Z,A,W,rho_diag,par,alpha,n_link,n_nonlink,type,method,gap_type)

    rho0p=par.rho0p;
    rho0n=par.rho0n;
    [logP_A,logP_S]=evalLikelihood(gap,n_link,n_nonlink,rho_diag,par,alpha,Z,type,method);
    noc=size(Z,1);
    J=size(A{1},1);    
    N=length(A);
    eN=ones(1,N);
    % step 1 select two observations i and j        
    ind1=ceil(J*rand);        
    ind2=ceil((J-1)*rand);
    if ind1<=ind2
       ind2=ind2+1;
    end
    clust1=find(Z(:,ind1));
    clust2=find(Z(:,ind2));

    if clust1==clust2 % Split   
        setS=find(sum(Z([clust1 clust2],:)));    
        setS=setdiff(setS,[ind1,ind2]);
        n_setS=length(setS);
        S_t=Z;
        S_t(clust1,:)=0;        
        comp=[clust1 noc+1];               
        S_t(comp(1),ind1)=1;
        S_t(comp(2),ind2)=1;  
        if strcmp(method,'Binary')
            rho_diag_t=[rho_diag; betarnd(rho0p(1)*eN,rho0n(1)*eN)];            
            logQ_trans_rho=sum(lnbetalike(rho_diag_t(end,:),rho0p(1)*eN,rho0n(1)*eN));
        else
            rho_diag_t=[rho_diag; gamrnd(rho0p(1)*eN,1/rho0n(1)*eN)];
            logQ_trans_rho=sum(lngammalike(rho_diag_t(end,:),rho0p(1)*eN,rho0n(1)*eN));
        end
        if strcmp(gap_type,'individual')
            gap_t=[gap; betarnd(eN,eN)];
        else
            gap_t=gap;
            % logQ_trans_gap=0;
        end
                
        % Reassign by restricted gibbs sampling                
        for rep=1:3
            [S_t,rho_diag_t,n_link_t,n_nonlink_t,gap_t,logQ_trans,comp]=Gibbs_sample_SBGC(gap_t,S_t,A,W,rho_diag_t,par,alpha,setS(randperm(n_setS)),type,method,gap_type,comp);                        
        end     
        [logP_A_t,logP_S_t]=evalLikelihood(gap_t,n_link_t,n_nonlink_t,rho_diag_t,par,alpha,S_t,type,method);               
        
        % Calculate Metropolis-Hastings ratio
        a_split=rand<exp(logP_A_t+logP_S_t-logP_A-logP_S-logQ_trans-logQ_trans_rho);         
        if a_split
           disp(['Splitting cluster ' num2str(clust1)])
           n_link=n_link_t;
           n_nonlink=n_nonlink_t;
           rho_diag=rho_diag_t;
           if strcmp(gap_type,'individual')
               gap=gap_t;
           end
           Z=S_t;
        end
    else % Merge                                     
        S_t=Z;
        S_t(clust1,:)=S_t(clust1,:)+S_t(clust2,:);
        setS=find(S_t(clust1,:));           
        S_t(clust2,:)=[];   
        rho_diag_t=rho_diag;
        rho_diag_t(clust2,:)=[];
        gap_t=gap;
        if strcmp(gap_type,'individual')            
            gap_t(clust2,:)=[];
        end
        if clust2<clust1
            clust1_t=clust1-1;
        else 
            clust1_t=clust1;
        end
        noc_t=noc-1;

        n_links=n_link;
        n_nonlinks=n_nonlink;
        Ap=rho0p(2)*ones(noc)+(rho0p(1)-rho0p(2))*eye(noc);
        An=rho0n(2)*ones(noc)+(rho0n(1)-rho0n(2))*eye(noc);
        for n=1:N
            n_links(:,:,n)=n_links(:,:,n)-Ap;
            n_nonlinks(:,:,n)=n_nonlinks(:,:,n)-An;
        end      
        n_link_t=n_links;
        n_link_t(clust1,:,:)=n_links(clust1,:,:)+n_links(clust2,:,:);
        n_link_t(:,clust1,:)=n_links(:,clust1,:)+n_links(:,clust2,:);        
        n_link_t(clust1,clust1,:)=n_link_t(clust1,clust1,:)+n_links(clust2,clust2,:);        
        n_link_t(clust2,:,:)=[];
        n_link_t(:,clust2,:)=[];
        n_nonlink_t=n_nonlinks;
        n_nonlink_t(clust1,:,:)=n_nonlinks(clust1,:,:)+n_nonlinks(clust2,:,:);
        n_nonlink_t(:,clust1,:)=n_nonlinks(:,clust1,:)+n_nonlinks(:,clust2,:);
        n_nonlink_t(clust1,clust1,:)=n_nonlink_t(clust1,clust1,:)+n_nonlinks(clust2,clust2,:);        
        n_nonlink_t(clust2,:,:)=[];
        n_nonlink_t(:,clust2,:)=[];
        Ap(clust2,:)=[];
        Ap(:,clust2)=[];
        An(clust2,:)=[];
        An(:,clust2)=[];
        for n=1:N
            n_link_t(:,:,n)=n_link_t(:,:,n)+Ap;
            n_nonlink_t(:,:,n)=n_nonlink_t(:,:,n)+An;
        end        
        
        [logP_A_t,logP_S_t]=evalLikelihood(gap_t,n_link_t,n_nonlink_t,rho_diag_t,par,alpha,S_t,type,method);                       

        % Split the merged cluster and calculate transition probabilties                
        % noc_tt=noc_t-1;
        setS=setdiff(setS,[ind1,ind2]);
        n_setS=length(setS);
        S_tt=S_t;
        S_tt(clust1_t,:)=0;        
        comp=[clust1_t noc_t+1];               
        S_tt(comp(1),ind1)=1;
        S_tt(comp(2),ind2)=1;                
        if strcmp(method,'Binary')            
            logQ_trans_rho=sum(lnbetalike(rho_diag(clust2,:),rho0p(1)*eN,rho0n(1)*eN));
        else            
            logQ_trans_rho=sum(lngammalike(rho_diag(clust2,:),rho0p(1)*eN,rho0n(1)*eN));
        end
        
        % Reassign by restricted gibbs sampling        
        for rep=1:2        
            [S_tt,rho_diag_tt,n_link_tt,n_nonlink_tt,gap_tt,logQ_trans,comp]=Gibbs_sample_SBGC(gap,S_tt,A,W,[rho_diag_t; rho_diag(clust2,:)],par,alpha,setS(randperm(n_setS)),type,method,gap_type,comp);               
        end
        Force=[1 2]*Z([clust1 clust2],:);        
        [S_tt,rho_diag_tt,n_link_tt,n_nonlink_tt,gap_tt,logQ_trans]=Gibbs_sample_SBGC(gap,S_tt,A,W,[rho_diag_t; rho_diag(clust2,:)],par,alpha,setS(randperm(n_setS)),type,method,gap_type,comp,Force);                        
        a_merge=rand<exp(logP_A_t+logP_S_t-logP_A-logP_S+logQ_trans+logQ_trans_rho);                 
        
        if a_merge
          disp(['Merging cluster ' num2str(clust1) ' with cluster ' num2str(clust2)])
          Z=S_t;          
          n_link=n_link_t;
          n_nonlink=n_nonlink_t;
          rho_diag=rho_diag_t;
          if strcmp(gap_type,'individual')
               gap=gap_t;
           end
        end
    end


% -------------------------------------------------------------------------  
function [logP_A,logP_S]=evalLikelihood(gap,n_link,n_nonlink,rho_diag,par,alpha,Z,type,method)
    
    rho0p=par.rho0p;
    rho0n=par.rho0n;
    N=size(n_link,3);
    [noc, J]=size(Z);    
    logP_A=0;            
    if strcmp(method,'Binary')
        clusterFun=@my_betaub;        
    else
        clusterFun=@my_gammaub;        
    end
    for n=1:N                
        if strcmp(type,'UnDirected')            
            if size(gap,2)>1
                logP_A=logP_A+sum(sum(triu(clusterFun(gap(:,n),n_link(:,:,n),n_nonlink(:,:,n),rho_diag(:,n),rho0p,rho0n))));
            else
                logP_A=logP_A+sum(sum(triu(clusterFun(gap,n_link(:,:,n),n_nonlink(:,:,n),rho_diag(:,n),rho0p,rho0n))));
            end
        else
            if size(gap,2)>0
                 logP_A=logP_A+sum(sum(clusterFun(gap(:,n),n_link(:,:,n),n_nonlink(:,:,n),rho_diag(:,n),rho0p,rho0n)));
            else
                logP_A=logP_A+sum(sum(clusterFun(gap,n_link(:,:,n),n_nonlink(:,:,n),rho_diag(:,n),rho0p,rho0n)));
            end
        end        
    end
    sumS=sum(Z,2);        
    logP_S=noc*log(alpha)+sum(gammaln(full(sumS)))-gammaln(J+alpha)+gammaln(alpha);    
    
 % -------------------------------------------------------------------------
 function [gap,logP]=MHsampleGap(gap,gap_type,n_link,n_nonlink,rho_diag,type,method,rho0p,rho0n,logP)
        
    if strcmp(method,'Binary')
        clusterFun=@my_betaub;        
    else
        clusterFun=@my_gammaub;        
    end    
    accept=0;    
    [noc,noct, N]=size(n_link);
    for rep=1:10  % sample for gap many times since this is inexpensive                                
        M=ceil(250*rand);
        if strcmp(gap_type,'same')
            if rem(rep,2)==1
                gap_t=betarnd(1,1);
                deltalogQ=0;
            elseif rep==2   
                a=M*gap;
                b=M-a;                                    
                gap_t=betarnd(a,b);                           
                a_new=M*gap_t;
                b_new=M-a_new;                                    
                deltalogQ=lnbetalike(gap_t,a,b)-lnbetalike(gap,a_new,b_new);                                        
            end                                
           Bnew=clusterFun(gap_t,n_link,n_nonlink,rho_diag,rho0p,rho0n);
           logPnew=0;
           for n=1:N
               if strcmp(type,'UnDirected')
                   logPnew=logPnew+sum(sum(triu(Bnew(:,:,n))));
               else
                   logPnew=logPnew+sum(sum(Bnew(:,:,n)));
               end
           end
            if rand<exp(logPnew-logP+deltalogQ) % notice log proposal = 0 hence ignored in MH move
               accept=accept+1;
               gap=gap_t;
               logP=logPnew;
            end                           
        else
            for n=1:N
                gap_t=gap(:,n);
                for d=randperm(noc)            
                   if strcmp(type,'UnDirected')
                       B=clusterFun(gap,n_link(:,d,n),n_nonlink(:,d,n),rho_diag(:,n),rho0p,rho0n,d);
                    else
                       B1=clusterFun(gap(:,n),n_link(:,d,n),n_nonlink(:,d,n),rho_diag(:,n),rho0p,rho0n,d);
                       B2=clusterFun(gap(:,n),n_link(d,:,n)',n_nonlink(d,:,n)',rho_diag(:,n),rho0p,rho0n,d);
                       B=B1+B2;
                       B(d)=B1(d);
                    end
                    if rem(rep,2)==1
                        gap_t(d)=betarnd(1,1);
                        deltalogQ=0;
                    else                        
                        a=M*gap(d,n);
                        b=M-a;                                    
                        gap_t(d)=betarnd(a,b);                           
                        a_new=M*gap_t(d);
                        b_new=M-a_new;                                    
                        deltalogQ=lnbetalike(gap_t(d),a,b)-lnbetalike(gap_t(d),a_new,b_new);                        
                    end          

                    if strcmp(type,'UnDirected')
                       Bnew=clusterFun(gap_t,n_link(:,d,n),n_nonlink(:,d,n),rho_diag(:,n),rho0p,rho0n,d);
                    else
                       B1=clusterFun(gap_t,n_link(:,d,n),n_nonlink(:,d,n),rho_diag(:,n),rho0p,rho0n,d);
                       B2=clusterFun(gap_t,n_link(d,:,n)',n_nonlink(d,:,n)',rho_diag(:,n),rho0p,rho0n,d);
                       Bnew=B1+B2;
                       Bnew(d)=B1(d);
                    end
                    if rand<exp(sum(Bnew)-sum(B)+deltalogQ) % notice log proposal = 0 hence ignored in MH move
                        accept=accept+1;
                        gap(d,n)=gap_t(d);
                    else
                        gap_t(d)=gap(d,n);
                    end                
                end 
            end
        end
    end
    gap(gap<eps)=eps; % Ensure numeric stability
    disp([' Accepted ' num2str(accept) ' samples for gap']);

    
% -------------------------------------------------------------------------
function [Z,rho_diag,n_link,n_nonlink,gap,logQ_trans,comp]=Gibbs_sample_SBGC(gap,Z,A,W,rho_diag,par,alpha,JJ,type,method,gap_type,comp,Force)        
    
    if nargin<13
        Force=[];
    end
    if nargin<12
        comp=[];
    end
    logQ_trans=0;

    rho0p=par.rho0p;
    rho0n=par.rho0n;
    switch method
        case 'Binary'
            clustFun=@my_betaub;
            diag_const=betaln(rho0p(1),rho0n(1));
        case 'Weighted'
            clustFun=@my_gammaub;
            diag_const=rho0p(1).*log(rho0n(1))-gammaln(rho0p(1));
    end
    T=3;
    N=length(A);
    [I,J]=size(A{1});
    if par.sameEta
        eN=par.N;    
    else
        eN=ones(1,N);    
    end
    t=0;   
    sumS=sum(Z,2);
    noc=length(sumS);    
            
    Ap=rho0p(1)*eye(noc)+rho0p(2)*(ones(noc)-eye(noc));
    An=rho0n(1)*eye(noc)+rho0n(2)*(ones(noc)-eye(noc));    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Count initial number of links and non-links between groups    
    if ~strcmp(type,'UnDirected')
        A2=cell(1,N);
        W2=cell(1,N);
    end
    SST=(sumS*sumS'-diag(sumS));    
    if par.sameEta
        SST=par.N*SST;        
    end
    rho_diag_t=cell(1,N);
    n_link=zeros(noc,noc,N);
    n_nonlink=zeros(noc,noc,N);
    for n=1:N        
        if strcmp(type,'UnDirected')
            switch method
                case 'Binary'   
                    SASt=Z*A{n}*Z';
                    SWSt=Z*W{n}*Z';
                    n_link(:,:,n)=SASt+SASt';    
                    n_link(:,:,n)=n_link(:,:,n)-0.5*diag(diag(n_link(:,:,n)));
                    n_link(:,:,n)=n_link(:,:,n)+Ap;
                    n_nonlink(:,:,n)=SST-SASt-SWSt-SASt'-SWSt';
                    n_nonlink(:,:,n)=n_nonlink(:,:,n)-0.5*diag(diag(n_nonlink(:,:,n)));
                    n_nonlink(:,:,n)=n_nonlink(:,:,n)+An;                        
                case 'Weighted'
                    SASt=Z*A{n}*Z';
                    SWSt=Z*W{n}*Z';
                    n_link(:,:,n)=SASt+SASt';    
                    n_link(:,:,n)=n_link(:,:,n)-0.5*diag(diag(n_link(:,:,n)));
                    n_link(:,:,n)=n_link(:,:,n)+Ap;
                    n_nonlink(:,:,n)=SST-SWSt-SWSt';
                    n_nonlink(:,:,n)=n_nonlink(:,:,n)-0.5*diag(diag(n_nonlink(:,:,n)));
                    n_nonlink(:,:,n)=n_nonlink(:,:,n)+An;                                            
            end  
            switch method
                case {'Binary'}
                    if par.sameEta
                        A{n}=A{n}+A{n}';
                    else
                        A{n}=logical(A{n}+A{n}');
                    end                    
                case {'Weighted'}
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
            switch method
                case 'Binary'                             
                    SASt=Z*A{n}*Z';
                    SWSt=Z*W{n}*Z';                    
                    n_link(:,:,n)=SASt+Ap;    
                    n_nonlink(:,:,n)=SST-SASt-SWSt+An;                                      
                case 'Weighted'
                    SASt=Z*A{n}*Z';
                    SWSt=Z*W{n}*Z';                    
                    n_link(:,:,n)=SASt+Ap;    
                    n_nonlink(:,:,n)=SST-SWSt+An;                                        
            end                
        end        
    end               
    cluster_eval=clustFun(gap,n_link,n_nonlink,rho_diag,rho0p,rho0n);       
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Main loop
    for k=JJ                                  
        t=t+1;
        if mod(t,5000)==0
            disp(['sampling ' num2str(t) ' out of ' num2str(J) ' nodes']);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Remove effect of s_k
        SA1k=zeros(noc,N);
        SW1k=zeros(noc,N);
        if ~strcmp('type','UnDirected')
            SA2k=zeros(noc,N);
            SW2k=zeros(noc,N);
        end
        for n=1:N
            SA1k(:,n)=Z*A{n}(:,k);
            SW1k(:,n)=Z*W{n}(:,k);            
            if ~strcmp(type,'UnDirected')
                SA2k(:,n)=Z*A2{n}(:,k);                
                SW2k(:,n)=Z*W2{n}(:,k);                                            
            end
        end          
        sumS=sumS-Z(:,k);   
        switch method
            case {'Binary'}
                nSA1k=sumS*eN-SA1k-SW1k;                
                if ~strcmp(type,'UnDirected')
                    nSA2k=sumS*eN-SA2k-SW2k;
                end
             case {'Weighted'}
                nSA1k=sumS*eN-SW1k;                
                if ~strcmp(type,'UnDirected')
                    nSA2k=sumS*eN-SW2k;
                end
        end
        d=find(Z(:,k));        
        % Remove link counts generated from assigment Z(:,k)
        if ~isempty(d)                       
            if strcmp(type,'UnDirected')
                n_link(:,d,:)=permute(n_link(:,d,:),[1 3 2])-SA1k;        
                if N==1
                    n_link(d,:)=n_link(d,:)-SA1k';
                else
                    n_link(d,:,:)=permute(n_link(d,:,:),[2 3 1])-SA1k;
                end
                n_nonlink(:,d,:)=permute(n_nonlink(:,d,:),[1 3 2])-nSA1k;               
                if N==1
                    n_nonlink(d,:)=n_nonlink(d,:)-nSA1k';                                      
                else
                    n_nonlink(d,:,:)=permute(n_nonlink(d,:,:),[2 3 1])-nSA1k;                                      
                end
                n_link(d,d,:)=permute(n_link(d,d,:),[3 1 2])+SA1k(d,:)';   
                n_nonlink(d,d,:)=permute(n_nonlink(d,d,:),[3 1 2])+nSA1k(d,:)';
            else
                n_link(:,d,:)=permute(n_link(:,d,:),[1 3 2])-SA1k;        
                if N==1
                    n_link(d,:)=n_link(d,:)-SA2k';
                else                    
                    n_link(d,:,:)=permute(n_link(d,:,:),[2 3 1])-SA2k;
                end
                n_nonlink(:,d,:)=permute(n_nonlink(:,d,:),[1 3 2])-nSA1k;               
                if N==1
                    n_nonlink(d,:)=n_nonlink(d,:)-nSA2k';                                                                     
                else
                    n_nonlink(d,:,:)=permute(n_nonlink(d,:,:),[2 3 1])-nSA2k;                                                                     
                end
            end                
        end
        Z(:,k)=0;               
        
        if isempty(comp) % Distinguish between restricted and non-restricted sampling
            % Non-restricted sampling
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Remove singleton cluster            
            if sumS(d)==0 
                v=1:noc;
                v(d)=[];
                length_d=length(d);
                d=[];
                noc=noc-length_d;                               
                P=sparse(1:noc,v,ones(1,noc),noc,noc+length_d);                        
                SA1k=P*SA1k;
                if ~strcmp(type,'UnDirected')
                    SA2k=P*SA2k; 
                end
                nSA1k=P*nSA1k; 
                if ~strcmp(type,'UnDirected')
                    nSA2k=P*nSA2k; 
                end
                Z=P*Z;                        
                sumS=sumS(v,1); 
                rho_diag=rho_diag(v,:);
                if strcmp(gap_type,'individual')
                    gap=gap(v,:);
                end
                n_link=n_link(v,v,:);            
                n_nonlink=n_nonlink(v,v,:);            
                cluster_eval=cluster_eval(v,v,:);                            
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Calculate probability for existing communties as well as proposal cluster
                        
            % Update cluster_eval without current node being assigned to any
            % clusters
            if ~isempty(d)
                cluster_eval(:,d,:)=clustFun(gap,n_link(:,d,:),n_nonlink(:,d,:),rho_diag,rho0p,rho0n,d); % removed the constant -betaln(Ap,An)))                               
                if strcmp(type,'UnDirected')
                     cluster_eval(d,:,:)=permute(cluster_eval(:,d,:),[2 1 3]);
                else
                    cluster_eval(d,:,:)=permute(clustFun(gap,permute(n_link(d,:,:),[2 1 3]),permute(n_nonlink(d,:,:),[2 1 3]),rho_diag,rho0p,rho0n,d),[2 1 3]); % removed the constant -betaln(Ap,An)))                                 
                end                                    
            end           
                        
            enoc=ones(1,noc);                            
            eT=ones(1,T);                            
            if strcmp(type,'UnDirected')                                
                % Update likelihood when assigning node to each of the
                % clusters
                link=zeros(noc,noc+T,N);
                link(:,1:noc,:)=n_link+permute(SA1k(:,:,enoc),[1 3 2]);
                link(:,noc+(1:T),:)=permute(SA1k(:,:,eT),[1 3 2])+rho0p(2);
                nolink=zeros(noc,noc+T,N);
                nolink(:,1:noc,:)=n_nonlink+permute(nSA1k(:,:,enoc),[1 3 2]);
                nolink(:,noc+(1:T),:)=permute(nSA1k(:,:,eT),[1 3 2])+rho0n(2);                                
                
                % Propose T new diagonal elements
                if strcmp(method,'Binary')
                    rho_diag_t=[rho_diag; betarnd(rho0p(1)*ones(T,N),rho0p(2)*ones(T,N))];                    
                else
                    rho_diag_t=[rho_diag; gamrnd(rho0p(1)*eT,1/rho0p(2)*eT)'];
                end
                if strcmp(gap_type,'individual')                    
                    gap_t=[gap; betarnd(ones(T,N),ones(T,N))];
                else
                    gap_t=gap;
                end                                
                
                % Evaluate likelihood without contribution of d^th cluster
                if N==1
                    %sum_cluster_eval_d=e*sum_cluster_eval-[sum(cluster_eval,1)'; zeros(T,n)];
                    sum_cluster_eval_d=-[sum(cluster_eval,1)'; zeros(T,n)];
                else
                    %sum_cluster_eval_d=e*sum_cluster_eval-[permute(sum(cluster_eval,1),[2 3 1]); zeros(T,n)];
                    sum_cluster_eval_d=-[permute(sum(cluster_eval,1),[2 3 1]); zeros(T,n)];
                end                
                
                % Evaluate likelihood with node assigned each cluster
                cluster_eval_d=clustFun(gap_t,link,nolink,rho_diag_t,rho0p,rho0n);  
                if N==1
                    sbeta=sum(cluster_eval_d,1)';
                else
                    sbeta=permute(sum(cluster_eval_d,1),[2 3 1]);
                end
                if strcmp(method,'Binary')
                    sbeta(noc+(1:T),:)=sbeta(noc+(1:T),:)...
                            +log(rho_diag_t(noc+(1:T),:))*(rho0p(1)-1)+log(1-rho_diag_t(noc+(1:T),:))*(rho0n(1)-1)-diag_const; % include within cluster likelihood     
                else
                    sbeta(noc+(1:T),:)=sbeta(noc+(1:T),:)...
                            +log(rho_diag_t(noc+(1:T),:))*(rho0p(1))-rho0n(1)-diag_const; % include within cluster likelihood
                end
                logQ=sum(sum_cluster_eval_d+sbeta,2); 
            else
                dbeta=zeros(noc,N);
                for n=1:N
                   dbeta(:,n)=diag(cluster_eval(:,:,n)); 
                end
                if N==1
                    % sum_cluster_eval_d=e*sum_cluster_eval-[sum(cluster_eval+cluster_eval',1)'-dbeta; zeros(T,n)];
                    sum_cluster_eval_d=-[sum(cluster_eval+cluster_eval',1)'-dbeta; zeros(T,n)];
                else
                    % sum_cluster_eval_d=e*sum_cluster_eval-[permute(sum(cluster_eval+permute(cluster_eval,[2 1 3]),1),[2 3 1])-dbeta; zeros(T,n)];
                    sum_cluster_eval_d=-[permute(sum(cluster_eval+permute(cluster_eval,[2 1 3]),1),[2 3 1])-dbeta; zeros(T,n)];
                end

                e_noc=ones(1,noc);
                link1=zeros(noc,noc+T,N);
                nolink1=zeros(noc,noc+T,N);
                link2=zeros(noc,noc+T,N);
                nolink2=zeros(noc,noc+1,N);                
                for n=1:N
                    link1(:,1:noc,n)=n_link(:,:,n)+SA1k(:,n*e_noc)+diag(SA2k(:,n));                             
                    nolink1(:,1:noc,n)=n_nonlink(:,:,n)+nSA1k(:,n*e_noc)+diag(nSA2k(:,n));
                    link2(:,1:noc,n)=n_link(:,:,n)'+SA2k(:,n*e_noc)+diag(SA1k(:,n));
                    nolink2(:,1:noc,n)=n_nonlink(:,:,n)'+nSA2k(:,n*e_noc)+diag(nSA1k(:,n));
                end
                link1(:,noc+(1:T),:)=permute(SA1k(:,:,eT),[1 3 2])+rho0p(2);
                nolink1(:,noc+(1:T),:)=permute(nSA1k(:,:,eT),[1 3 2])+rho0n(2);                                                
                link2(:,noc+(1:T),:)=permute(SA2k(:,:,eT),[1 3 2])+rho0p(2);                                
                nolink2(:,noc+(1:T),:)=permute(nSA2k(:,:,eT),[1 3 2])+rho0n(2);                                                                
                
                % Propose T new diagonal elements
                if strcmp(method,'Binary')
                    rho_diag_t=[rho_diag; betarnd(rho0p(1)*ones(T,N),rho0p(2)*ones(T,N))];
                else
                    rho_diag_t=[rho_diag; gamrnd(rho0p(1)*eT,1/rho0p(2)*eT)'];
                end
                if strcmp(gap_type,'individual')                    
                    gap_t=[gap; betarnd(ones(T,N),ones(T,N))];
                else
                    gap_t=gap;
                end
                
                cluster_eval_d=clustFun(gap_t,link1,nolink1,rho_diag_t,rho0p,rho0n);                     
                cluster_eval_d_2=clustFun(gap_t,link2,nolink2,rho_diag_t,rho0p,rho0n);

                dbeta=zeros(noc+T,N);
                for n=1:N
                   dbeta(1:noc,n)=diag(cluster_eval_d(:,1:noc,n)); 
                end
                if N==1
                    sbeta=sum(cluster_eval_d+cluster_eval_d_2,1)'-dbeta;
                else
                    sbeta=permute(sum(cluster_eval_d+cluster_eval_d_2,1),[2 3 1])-dbeta;
                end
                %sbeta(noc+(1:T),:)=sbeta(noc+(1:T),:)...                    
                %    +log(rho_diag_t(noc+(1:T),:))*(rho0p(1)-1)+log(1-rho_diag_t(noc+(1:T),:))*(rho0n(1)-1); % include within cluster likelihood
                logQ=sum(sum_cluster_eval_d+sbeta,2); 
            end 
                
                    
            % Sample from posterior                  
            QQ=exp(logQ-max(logQ));
            weight=[sumS; alpha/T*eT'];
            QQ=weight.*QQ;                            
            ind=find(rand<full(cumsum(QQ/sum(QQ))),1,'first');                                         
            if ind>noc                                    
                noc=noc+1;
                Z(noc,k)=1;   
                rho_diag=[rho_diag; rho_diag_t(ind,:)];                
                if strcmp(gap_type,'individual')                    
                    gap=[gap; gap_t(ind,:)];
                end                        
                sumS(noc,1)=0;
                n_link(:,noc,:)=rho0p(2);
                n_link(noc,:,:)=rho0p(2);            
                n_link(noc,noc,:)=rho0p(1);            
                n_nonlink(:,noc,:)=rho0n(2);                     
                n_nonlink(noc,:,:)=rho0n(2);                     
                n_nonlink(noc,noc,:)=rho0n(1);                                 
                cluster_eval(:,noc,:)=0;    
                cluster_eval(noc,:,:)=0;                
                cluster_eval_d1=permute(cluster_eval_d(:,ind,:),[1 3 2]);
                cluster_eval_d1(noc,:)=0;                                
                SA1k(noc,:)=0;
                nSA1k(noc,:)=0;              
                if ~strcmp(type,'UnDirected')
                    if N==1
                        cluster_eval_d2=cluster_eval_d_2(:,ind);
                    else
                        cluster_eval_d2=permute(cluster_eval_d_2(:,ind,:),[1 3 2]);                                                
                    end                    
                    cluster_eval_d2(noc,:)=0;
                    SA2k(noc,:)=0;
                    nSA2k(noc,:)=0;          
                end                   
                ind=noc;
            else
                Z(ind,k)=1;   
                cluster_eval_d1=permute(cluster_eval_d(1:noc,ind,:),[1 3 2]);
                if ~strcmp(type,'UnDirected')
                    cluster_eval_d2=permute(cluster_eval_d_2(1:noc,ind,:),[1 3 2]);
                end                
            end                        
        else            
            % Calculate probability for existing communties as well as proposal cluster                                                            
            if ~isempty(d)
                cluster_eval(:,d,:)=clustFun(gap,n_link(:,d,:),n_nonlink(:,d,:),rho_diag,rho0p,rho0n,d); % removed the constant -betaln(Ap,An)))                               
                if strcmp(type','UnDirected')
                    cluster_eval(d,:,:)=squeeze(cluster_eval(:,d,:));
                else
                    cluster_eval(d,:,:)=permute(clustFun(gap,permute(n_link(d,:,:),[2 1 3]),permute(n_nonlink(d,:,:),[2 1 3]),rho_diag,rho0p,rho0n,d),[2 1 3]); % removed the constant -betaln(Ap,An)))                               
                end
            end
            %for n=1:N
            %    if strcmp(type,'UnDirected')
            %        sum_cluster_eval(n)=sum(sum(triu(cluster_eval(:,:,n))));                        
            %    else
            %        sum_cluster_eval(n)=sum(sum(cluster_eval(:,:,n)));                        
            %    end
            %end
            e=ones(2,1);
            if strcmp(type,'UnDirected')
                if N==1
                  %  sum_cluster_eval_d=e*sum_cluster_eval-sum(cluster_eval(:,comp))';                        
                    sum_cluster_eval_d=-sum(cluster_eval(:,comp))';                        
                else
                  %  sum_cluster_eval_d=e*sum_cluster_eval-permute(sum(cluster_eval(:,comp,:)),[2 3 1]);                        
                    sum_cluster_eval_d=-permute(sum(cluster_eval(:,comp,:)),[2 3 1]);                        
                end
                link=n_link(:,comp,:)+permute(SA1k(:,:,e),[1 3 2]);                        
                nolink=n_nonlink(:,comp,:)+permute(nSA1k(:,:,e),[1 3 2]);            
                cluster_eval_d1=clustFun(gap,link,nolink,rho_diag,rho0p,rho0n,comp);
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
                    % sum_cluster_eval_d=e*sum_cluster_eval-(sum(cluster_eval(:,comp)+cluster_eval(comp,:)')'-dbeta);                        
                    sum_cluster_eval_d=-(sum(cluster_eval(:,comp)+cluster_eval(comp,:)')'-dbeta);                        
                else
                    % sum_cluster_eval_d=e*sum_cluster_eval-(permute(sum(cluster_eval(:,comp,:)+permute(cluster_eval(comp,:,:),[2 1 3])),[2 3 1])-dbeta);                        
                     sum_cluster_eval_d=-(permute(sum(cluster_eval(:,comp,:)+permute(cluster_eval(comp,:,:),[2 1 3])),[2 3 1])-dbeta);                        
                end
                e_noc=ones(1,2);
                link1=n_link;
                nolink1=n_nonlink;
                link2=permute(n_link,[2 1 3]);
                nolink2=permute(n_nonlink,[2 1 3]);
                for n=1:N
                    link1(:,comp,n)=link1(:,comp,n)+SA1k(:,n*e_noc);                                        
                    link1(comp,comp,n)=link1(comp,comp,n)+diag(SA2k(comp,n));               
                    nolink1(:,comp,n)=nolink1(:,comp,n)+nSA1k(:,n*e_noc);                            
                    nolink1(comp,comp,n)=nolink1(comp,comp,n)+diag(nSA2k(comp));                                
                    link2(:,comp,n)=link2(:,comp,n)+SA2k(:,n*e_noc);                                        
                    link2(comp,comp,n)=link2(comp,comp,n)+diag(SA1k(comp,n));
                    nolink2(:,comp,n)=nolink2(:,comp,n)+nSA2k(:,n*e_noc);                            
                    nolink2(comp,comp,n)=nolink2(comp,comp,n)+diag(nSA1k(comp,n));
                end                
                cluster_eval_d1=clustFun(gap,link1(:,comp,:),nolink1(:,comp,:),rho_diag,rho0p,rho0n,comp);     
                cluster_eval_d2=clustFun(gap,link2(:,comp,:),nolink2(:,comp,:),rho_diag,rho0p,rho0n,comp);     
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

            % Sample from posterior  
            QQ=exp(logQ-max(logQ));
            weight=sumS(comp);
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
            if strcmp(type,'UnDirected')
                cluster_eval_d1=cluster_eval_d1(:,ind,:);                
            else
                cluster_eval_d1=permute(cluster_eval_d1(:,ind,:),[1 3 2]);
                cluster_eval_d2=permute(cluster_eval_d2(:,ind,:),[1 3 2]);
            end            
            ind=comp(ind);            
        end
                        
        % Re-enter effect of new s_k        
        sumS=sumS+Z(:,k);

        if strcmp(type,'UnDirected')
            n_link(:,ind,:)=permute(n_link(:,ind,:),[1 3 2])+SA1k;
            if N==1
                n_link(ind,:)=n_link(ind,:)+SA1k';
            else
                n_link(ind,:,:)=permute(n_link(ind,:,:),[2 3 1])+SA1k;
            end
            n_link(ind,ind,:)=permute(n_link(ind,ind,:),[3 1 2])-SA1k(ind,:)';
            n_nonlink(:,ind,:)=permute(n_nonlink(:,ind,:),[1 3 2])+nSA1k;        
            if N==1
                n_nonlink(ind,:)=n_nonlink(ind,:)+nSA1k';                
            else
                n_nonlink(ind,:,:)=permute(n_nonlink(ind,:,:),[2 3 1])+nSA1k;                
            end
            n_nonlink(ind,ind,:)=permute(n_nonlink(ind,ind,:),[3 1 2])-nSA1k(ind,:)';                
            cluster_eval(:,ind,:)=cluster_eval_d1;
            cluster_eval(ind,:,:)=cluster_eval_d1;
        else
            n_link(:,ind,:)=permute(n_link(:,ind,:),[1 3 2])+SA1k;                
            if N==1
                n_link(ind,:)=n_link(ind,:)+SA2k';            
            else
                n_link(ind,:,:)=permute(n_link(ind,:,:),[2 3 1])+SA2k;            
            end
            n_nonlink(:,ind,:)=permute(n_nonlink(:,ind,:),[1 3 2])+nSA1k;        
            if N==1
                n_nonlink(ind,:)=n_nonlink(ind,:)+nSA2k';                            
            else
                n_nonlink(ind,:,:)=permute(n_nonlink(ind,:,:),[2 3 1])+nSA2k;                            
            end
            cluster_eval(:,ind,:)=cluster_eval_d1;
            if N==1
                cluster_eval(ind,:)=cluster_eval_d2';
            else
                cluster_eval(ind,:,:)=cluster_eval_d2;
            end
        end                                                    
                
        % Remove empty clusters        
        if ~all(sumS)
            d=find(sumS==0);
            ind_d=find(d<comp);
            comp(ind_d)=comp(ind_d)-1;
            v=1:noc;
            v(d)=[];
            rho_diag=rho_diag(v,:);
            if strcmp(gap_type,'individual')
                gap=gap(v,:);
            end
            noc=noc-1;                               
            P=sparse(1:noc,v,ones(1,noc),noc,noc+1);                        
            Z=P*Z;                        
            sumS=sumS(v,1);    
            n_link=n_link(v,v,:);
            n_nonlink=n_nonlink(v,v,:);
            cluster_eval=cluster_eval(v,v,:);                            
        end           
        
    end                          
    
% --------------------------------------
 function B=my_betaub(gap,n_link,n_nonlink,rho_diag,rho0p,rho0n,dd)
        % n_link        D x DD x N matrix
        % n_nonlink     D x DD x N matrix
        % rho_diag      D x N matrix 
        % dd             indices corresponding to DD
                
        [D,DD,N]=size(n_link);                             
        beta_diag_const=betaln(rho0p(1),rho0n(1));
        
        if nargin<7
            if numel(gap)==1
                rho_min=gap*minComDens(rho_diag);           
                E=ones(size(rho_diag));            
                BetaIncp=betainc(gap*rho_diag,rho0p(2)*E,rho0n(2)*E);
                BetaIncp(BetaIncp==0)=1e-323;
                prior_diag=betaln(rho0p(2),rho0p(2))+log(BetaIncp);    
            else
                rho_min=minComDens(gap.*rho_diag);           
                E=ones(size(rho_diag));           
                BetaIncp=betainc(gap.*rho_diag,rho0p(2)*E,rho0n(2)*E);
                BetaIncp(BetaIncp==0)=1e-323;
                prior_diag=betaln(rho0p(2),rho0p(2))+log(BetaIncp);    
            end
            
            prior_min=minComDens(prior_diag);
            BetaInc=betainc(rho_min(1:size(n_link,1),:,:),n_link,n_nonlink);
            BetaInc(BetaInc==0)=1e-323;
            logBetaInc=log(BetaInc);
            B = betaln(n_link,n_nonlink)+logBetaInc-prior_min(1:size(n_link,1),:,:);        
            for n=1:N
                if D==DD
                    dn_link=diag(n_link(:,:,n));
                    dn_nonlink=diag(n_nonlink(:,:,n));
                    ii=find(diag(B(:,:,n))==-Inf);
                    for t=1:length(ii)
                        B(ii(t),ii(t))=0;
                    end
                    B(:,:,n)=B(:,:,n)-diag(diag(B(:,:,n)))+diag(log(rho_diag(:,n)).*(dn_link-1)+log(1-rho_diag(:,n)).*(dn_nonlink-1)-beta_diag_const);           
                else
                    dn_link=diag(n_link(1:D,1:D,n));
                    dn_nonlink=diag(n_nonlink(1:D,1:D,n));
                    ii=find(diag(B(1:D,1:D,n))==-Inf);
                    for t=1:length(ii)
                        B(ii(t),ii(t))=0;
                    end
                    B(1:D,1:D,n)=B(1:D,1:D,n)-diag(diag(B(1:D,1:D,n)))+diag(log(rho_diag(1:D,n)).*(dn_link-1)+log(1-rho_diag(1:D,n)).*(dn_nonlink-1)-beta_diag_const);           
                end
            end                   
        else     
            if numel(gap)==1
                rho_min=gap*minComDens(rho_diag,dd);   
                E=ones(size(rho_diag));
                BetaIncp=betainc(gap*rho_diag,rho0p(2)*E,rho0n(2)*E);
                BetaIncp(BetaIncp==0)=1e-323;
                prior_diag=betaln(rho0p(2),rho0p(2))+log(BetaIncp);    
            else
                rho_min=minComDens(gap.*rho_diag,dd);          
                E=ones(size(rho_diag));
                BetaIncp=betainc(gap.*rho_diag,rho0p(2)*E,rho0n(2)*E);
                BetaIncp(BetaIncp==0)=1e-323;                
                prior_diag=betaln(rho0p(2),rho0p(2))+log(BetaIncp);    
            end
            
            prior_min=minComDens(prior_diag,dd);   
            BetaInc=betainc(rho_min(1:size(n_link,1),:,:),n_link,n_nonlink);
            BetaInc(BetaInc==0)=1e-323;
            logBetaInc=log(BetaInc);
            B = betaln(n_link,n_nonlink)+logBetaInc-prior_min(1:size(n_link,1),:,:);                    
            if dd<=D                            
                for n=1:N
                    dn_link=diag(n_link(dd,:,n));
                    dn_nonlink=diag(n_nonlink(dd,:,n));
                    ii=find(diag(B(dd,:,n))==-Inf);
                    for t=1:length(ii)
                        B(dd(ii(t)),ii(t),n)=0;
                    end
                    B(dd,:,n)=B(dd,:,n)-diag(diag(B(dd,:,n)))+diag(log(rho_diag(dd,n)).*(dn_link-1)+log(1-rho_diag(dd,n)).*(dn_nonlink-1)-beta_diag_const);           
                end
            end            
        end
        
% --------------------------------------
 function B=my_gammaub(gap,n_link,n_nonlink,rho_diag,rho0p,rho0n,dd)
        % n_link        D x DD x N matrix
        % n_nonlink     D x DD x N matrix
        % rho_diag      D x N matrix 
        % dd             indices corresponding to DD
        
        [D,DD,N]=size(n_link);             
        gamma_diag_const=rho0p(1).*log(rho0n(1))-gammaln(rho0p(1));
        if nargin<7
            if numel(gap)==1
                rho_min=gap*minComDens(rho_diag);                    
                E=ones(size(rho_diag));
                GamIncp=gammainc(gap*rho_diag*rho0n(2),rho0p(2)*E);
                GamIncp(GamIncp==0)=1e-323;
                prior_diag=-rho0p(2).*log(rho0n(2))+gammaln(rho0p(2))+log(GamIncp);                
            else
                rho_min=minComDens(gap.*rho_diag);                    
                E=ones(size(rho_diag));
                GamIncp=gammainc(gap.*rho_diag*rho0n(2),rho0p(2)*E);
                GamIncp(GamIncp==0)=1e-323;
                prior_diag=-rho0p(2).*log(rho0n(2))+gammaln(rho0p(2))+log(GamIncp);                
            end            
            prior_min=minComDens(prior_diag);            
            GamInc=gammainc(rho_min(1:size(n_link,1),:,:).*n_nonlink,n_link);
            GamInc(GamInc==0)=1e-323;
            logGamInc=log(GamInc);            
            B = -n_link.*log(n_nonlink)+gammaln(n_link)+logGamInc-prior_min(1:size(n_link,1),:,:);        
            for n=1:N
                if D==DD
                    dn_link=diag(n_link(:,:,n));
                    dn_nonlink=diag(n_nonlink(:,:,n));
                    ii=find(diag(B(:,:,n))==-Inf);
                    for t=1:length(ii)
                        B(ii(t),ii(t),n)=0;
                    end
                    B(:,:,n)=B(:,:,n)-diag(diag(B(:,:,n)))+diag(log(rho_diag(:,n)).*(dn_link-1)-rho_diag(:,n).*dn_nonlink-gamma_diag_const);           
                else
                    dn_link=diag(n_link(1:D,1:D,n));
                    dn_nonlink=diag(n_nonlink(1:D,1:D,n));
                    ii=find(diag(B(1:D,1:D,n))==-Inf);
                    for t=1:length(ii)
                        B(ii(t),ii(t))=0;
                    end
                    B(1:D,1:D,n)=B(1:D,1:D,n)-diag(diag(B(1:D,1:D,n)))+diag(log(rho_diag(1:D,n)).*(dn_link-1)-rho_diag(1:D,n).*dn_nonlink-gamma_diag_const);           
                end
            end                   
        else
            if numel(gap)==1
                rho_min=gap*minComDens(rho_diag,dd);                    
                E=ones(size(rho_diag));
                GamIncp=gammainc(gap*rho_diag*rho0n(2),rho0p(2)*E);
                GamIncp(GamIncp==0)=1e-323;
                prior_diag=-rho0p(2).*log(rho0n(2))+gammaln(rho0p(2))+log(GamIncp);                
            else
                rho_min=minComDens(gap.*rho_diag,dd);                    
                E=ones(size(rho_diag));                
                GamIncp=gammainc(gap.*rho_diag*rho0n(2),rho0p(2)*E);
                GamIncp(GamIncp==0)=1e-323;
                prior_diag=-rho0p(2).*log(rho0n(2))+gammaln(rho0p(2))+log(GamIncp);                
            end            
            
            prior_min=minComDens(prior_diag,dd);     
            GamInc=gammainc(rho_min(1:size(n_link,1),:,:).*n_nonlink,n_link);
            GamInc(GamInc==0)=1e-323;
            logGamInc=log(GamInc);
            B = -n_link.*log(n_nonlink)+gammaln(n_link)+logGamInc-prior_min(1:size(n_link,1),:,:);
            if dd<=D  
                for n=1:N
                    dn_link=diag(n_link(dd,:,n));
                    dn_nonlink=diag(n_nonlink(dd,:,n));
                    ii=find(diag(B(dd,:,n))==-Inf);
                    for t=1:length(ii)
                        B(dd(ii(t)),ii(t),n)=0;
                    end
                    B(dd,:,n)=B(dd,:,n)-diag(diag(B(dd,:,n)))+diag(log(rho_diag(dd,n)).*(dn_link-1)-rho_diag(dd,n).*dn_nonlink-gamma_diag_const);           
                end
            end
        end

        
    

        
% ----------------------------------------------------
function rho_min=minComDens(rho_diag,d)
[noc,N]=size(rho_diag);
if nargin<2
    ed=ones(1,noc);
    rho_cc=zeros(noc,noc,N,2);
else
    ed=ones(1,length(d));
    enoc=ones(noc,1);
    rho_cc=zeros(noc,length(d),N,2);
end
rho_c=permute(rho_diag(:,:,ed),[1 3 2]);
rho_cc(:,:,:,1)=rho_c;
if nargin<2
    rho_cc(:,:,:,2)=permute(rho_c,[2 1 3]);   
else
    rho_cc(:,:,:,2)=permute(rho_diag(d,:,enoc),[3 1 2]);
end
rho_min=min(rho_cc,[],4);


% ----------------------------------------------------
function logL=lnbetalike(x,a,b)    
    logL=gammaln(a+b)-gammaln(a)-gammaln(b)+(a-1).*log(x+1e-323)+(b-1).*log(1-x+1e-323);

% ----------------------------------------------------
function logL=lngammalike(x,a,b)    
    logL=a.*log(b)-gammaln(a)+(a-1).*log(x)-b.*x;
            
%-----------------------------------------------------
    function noise_est=my_noise_est_beta(gap,n_link,n_nonlink)
        noise_est=nan(size(n_link));
        steps=10000;
        for t=1:length(n_link)
            x=linspace(1e-320,gap(t),steps);
            q=(n_nonlink(t)-1)*log(1-x);
            logP1=(n_link(t)+1-1)*log(x)+q;
            logP2=(n_link(t)-1).*log(x)+q;        
            maxlogP1=max(logP1);
            logP1=logP1-maxlogP1;
            logP2=logP2-maxlogP1;
            noise_est(t)=trapz(exp(logP1))/trapz(exp(logP2));
        end
        
%-----------------------------------------------------
    function noise_est=my_noise_est_gamma(gap,n_link,n_nonlink)
        noise_est=nan(size(n_link));
        steps=10000;
        for t=1:length(n_link)
            x=linspace(1e-320,gap(t),steps);        
            logP1=(n_link(t)+1-1)*log(x)-n_nonlink(t)*x;
            logP2=(n_link(t)-1).*log(x)-n_nonlink(t)*x;        
            maxlogP1=max(logP1);
            logP1=logP1-maxlogP1;
            logP2=logP2-maxlogP1;
            noise_est(t)=trapz(exp(logP1))/trapz(exp(logP2));
        end
        
        