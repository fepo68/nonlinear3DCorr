%% Function for Gradient Based Optimization of the Projection Matrix UCM-LVM

% H.F. Garcia 2016

function [f,g] = projectionMatrixOptFCN(x,params,X,S,d)
% Data setup
% [m,n] = size(Data.Wd);
Saux = S{d};
Xd = X{d};
Wd = reshape(x,params.Md(d),params.K);
k =rank(Wd);


% Function Value 
gammaVal = params.gamma;
J = params.J;
K = params.K;

auxNj = 0;

for j = 1:J
    Nj = Ndot_j(S,j,params);
    for i = 1:Nj-1
    auxNj = log(i)+auxNj;
    end
end
auxN = 0;
for n = 0:params.N-1
    auxN = auxN + log(gammaVal+n);
end
    

auxSpart = J*log(gammaVal)+auxNj-auxN;

auxXpart = - params.auxSumD*log(2*pi)+((K*J)/2)*log(params.r)+params.a*log(params.b)...
    -params.ap*log(params.bp)+(gammaln(params.ap))-(gammaln(params.a));

auxdetCj  = 0;

for j = 1:J
    auxdetCj = auxdetCj +0.5*log(det(inv(params.invCj(:,:,j))));
end

f = auxSpart+auxXpart+auxdetCj;
% First derivative computed in matrix form for log p(X,S|W,a,b,r,gamma)
% with respect to Wd

auxCj = zeros(K);
auxCjpart2 = zeros(params.Md(d),K);
for j = 1:J
    % Loop for Cj
    Ndj = Ndjval(S,j,d);
    auxCj = auxCj + Ndj*inv(params.invCj(:,:,j));
    % Loop for Wd mu_j and xdn
    % Fin the number of objects in the domain d assigned to the cluster j
    [posObj,nObj_j] = find(Saux == j);
    % Rigth side of the equation 18
    sumXdn_j = sum(Xd(posObj,:),1)';
    auxCjpart2 = auxCjpart2 + Ndj*Wd*(params.mu_j(:,j)*params.mu_j(:,j)')-(sumXdn_j*params.mu_j(:,j)');
end
G = -Wd*auxCj-(params.ap/params.bp)*auxCjpart2;
g = G(:);

