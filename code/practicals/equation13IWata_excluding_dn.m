%% Second factor parameters depicted in (4) but excluding dn object
% see bottom factor of the equation 13
function [ap_no_dn,bp,mu_j,invCj] = equation13IWata_excluding_dn(X,S,W,params,de,ne)
a = params.a;
b= params.b;
r = params.r;
gamma = params.gamma;
K = params.K;
J = params.J;
D = params.D;
auxN = params.auxN;

nNd = params.Nd;
nNd(de) = nNd(de)-1; 


ap_no_dn = a + sum(nNd.*params.Md)/2;
% ap_no_dn = a + (auxN-1)/2;
invCj = zeros(K,K,J);
mu_j = zeros(K,J);
for j = 1:J
    % For covariances C_{j_\dn}^-1
    auxCj = zeros(K);
    for d = 1:D
        auxS = S{d};
        [~,n2]=find(auxS ==j);
        if de == d
            if S{de}(ne) == j
            Ndj = length(n2)-1;
            else
                Ndj = length(n2);
            end
        else
            Ndj = length(n2);
        end
        Wd = W{d};
        auxCj = auxCj + Ndj*(Wd'*Wd);
    end
    invCj(:,:,j) = auxCj+r*eye(K);
    % For means <mu>
    aux_mj = zeros(K,1);
    for d = 1:D
        auxS = S{d};
        if de ~= d
            [n1,n2]=find(auxS ==j); % Find objects that
            % has assigned the j-th cluster
        else
            [n1,n2]=find(auxS ==j); % Find objects that
            % has assigned the j-th cluster, except the dn object
            [i_ne,~] = find(n1==ne);
            n1(i_ne) = [];
            n2(i_ne) = [];
        end
        
        Wd = W{d};
        Xd = [X{d}]'; % Because objects are Md*Nd
        aux_xdn =sum(Xd(:,n1),2);
        aux_mj = aux_mj + Wd'*aux_xdn;
    end
    mu_j(:,j) = invCj(:,:,j)\aux_mj;
    
end

% For b' eq(6) but excluding the dn object
aux_xdn = 0;
for d = 1:D
    auxS = S{d};
    if d ==de
%         [i_ne,~] = find(auxS==ne);
        Xd = [X{d}]';
        Xd(:,ne) = [];
        [Md,Nd] = size(Xd);
    else
        Xd = [X{d}]';
        [Md,Nd] = size(Xd);
    end
    
    aux = sum(diag(Xd'*Xd));
    aux_xdn = aux_xdn + aux;
%     for n = 1:Nd
%         %         if ne ~= n
%         xdn = Xd(:,n);
%         aux_xdn = aux_xdn +(xdn'*xdn);
%         %         end
%     end
end
aux_muCj = 0;
for j = 1:J
    aux_muCj = aux_muCj + (mu_j(:,j)'*invCj(:,:,j)*mu_j(:,j));
end
bp = b + (1/2)*aux_xdn-(1/2)*aux_muCj;