function [Sclust,Y,Fd,X] = createToyMultiviewWMMCluster(Nd,Dd,V,K,Q,beta)

% Draw samples from the Multiview Warped Mixture Model
% Inputs:
% Nd -> array where each element contains the number of objects per view
% Dd -> array where each element contains the number of outputs per view
% V -> Number of views
% K -> Number of Clusters
% Q -> Dimensionality of the Latent Space

GEM = false;
% 1. Draw mixture weights
if GEM == true
    % Initialize lambda with a Stic Breaking Process
else
    % All mixture weigths have the same probability
    lambda = ones(1,K)/K;
end
% 2. For each component c = 1,..., infty
R = zeros(Q,Q,K);
mu = zeros(Q,K);
S = eye(Q);
nu = 15; % degree of freedom for the Wishart dist
u = rand(Q,K);
r = 0.1;
for k = 1:K
    % Draw preciosion Rc
    R(:,:,k) = wishrnd(S,nu);
    % Draw mean
    mu(:,k) = mvnrnd(u(:,k),(r*R(:,:,k))\eye(Q));
    
end
Z = [];
X = {};
Sclust = {};
for v = 1:V
    Xd_aux = zeros(Nd(v),Q);
    Zv = mnrnd(1,lambda,Nd(v));
    [~, Sclust{v}] = max(Zv');
    for n = 1:Nd(v)
        [~,idZn] = max(Zv(n,:));
        Xd_aux(n,:) = mvnrnd(mu(:,idZn),inv(R(:,:,idZn)));
    end
    X{v} = Xd_aux;
    %     [~,idZn] = max(Zv');
    %     idZn = idZn';
    %     Xd_aux = mvnrnd(mu(:,idZn),inv(R(:,:,idZn)),Nd(v));
    
    %% Draw functions
    if v == 1
        kern = kernCreate(Xd_aux, 'rbf');
        kern.inverseWidth = 0.1^2;
        kern.variance = 1;
    else
        kern = kernCreate(Xd_aux, 'rbfard');
        kern.inverseWidth = 0.01;
        kern.variance = 0.8;
    end
    K = kernCompute(kern, Xd_aux);
%     imagesc(K);
    % need to take the real part of the sample as the kernel is numerically less than full rank
    % Sample fdm using a fdm(z)~GP(0,K)
    Fd{v} = real(gsamp(zeros(1, size(Xd_aux, 1)), K, Dd(v)))';
    for d = 1:Dd
    end
end

% Draw Observations
Y = {};
for v = 1:V
    Fdaux = Fd{v};
    Yaux = zeros(Nd(v),Dd(v));
    for n = 1:Nd(v);
        Yaux(n,:) = mvnrnd(Fdaux(n,:),(beta^(-1))*eye(Dd(v)));
    end
    Y{v} = Yaux;
end

