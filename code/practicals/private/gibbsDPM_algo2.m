function c_st = gibbsDPM_algo2(y, hyperG0, alpha, niter, doPlot)

% Gibbs sampler for Dirichlet Process Mixtures of Gaussians
% The base distribution G0 is Normal inverse Wishart of parameters given by
% hyperG0
% Reference: Algorithm 2 of Neal

if doPlot
    figure('name','Gibbs sampling for DPM');
    colormap('default')
    cmap = colormap;
end


[p, n] = size(y);
c_st = zeros(n, niter/2);
U_mu = zeros(p, n);
U_Sigma = zeros(p, p, n);

% U_SS is a structure array where U_SS(k) contains the sufficient
% statistics associated to cluster k
U_SS = struct('mu', cell(n, 1), 'kappa', cell(n, 1), ...
    'nu', cell(n, 1), 'lambda', cell(n, 1));

m=zeros(1,200);
c = zeros(n, 1);
% Initialisation
for k=1:n
    c(k) = ceil(30*rand); % Sample new allocation uniform
    m(c(k)) = m(c(k)) + 1;
    if m(c(k))>1
        U_SS(c(k)) = update_SS(y(:,k), U_SS(c(k)));
    else
        U_SS(c(k)) = update_SS(y(:,k), hyperG0);
    end
end
ind = unique(c);
for j=1:length(ind)
    [U_mu(:, ind(j)), U_Sigma(:, :, ind(j))] = normalinvwishrnd(U_SS(ind(j)));
end


% Iterations
for i=2:niter
    % Update cluster assignments c
    for k=1:n        
        m(c(k)) = m(c(k)) - 1;
        U_SS(c(k)) = downdate_SS(y(:,k),U_SS(c(k)));
        c(k) = sample_c(m, alpha, y(:,k), hyperG0, U_mu, U_Sigma);
        m(c(k)) = m(c(k)) + 1;
        if m(c(k))>1
            U_SS(c(k)) = update_SS(y(:,k), U_SS(c(k)));
        else
            U_SS(c(k)) = update_SS(y(:,k), hyperG0);
           [U_mu(:, c(k)), U_Sigma(:, :, c(k))] = normalinvwishrnd(U_SS(c(k)));
        end       
      
        if doPlot==1
            some_plot(y, U_mu, m, c, k, i, cmap)
        end
    end
    % Update cluster locations U
    ind = find(m);
    for j=1:length(ind)
        [U_mu(:, ind(j)), U_Sigma(:, :, ind(j))] = normalinvwishrnd(U_SS(ind(j)));
    end
    
    fprintf('Iteration %d/%d\n', i, niter)
    fprintf('%d clusters\n\n', length(unique(c)))
    
    if doPlot==2
        some_plot(y, U_mu, m, c, k, i, cmap)
    end

    if i>niter/2
        c_st(:, i-niter/2) = c;
    end

end
end


function K = sample_c(m, alpha, z, hyperG0, U_mu, U_Sigma)

    c = find(m~=0); % gives indices of non-empty clusters
    r = sum(m);
    n = m(c).*exp(loggausspdf(repmat(z, 1, length(c))', U_mu(:, c)', U_Sigma(:, :, c))');

    n0 = pred(z, hyperG0);
    const = sum(n) + alpha*n0;

    p0 = alpha*n0/const; % probability of sampling a new item

    u=rand(1);
    if u<p0
        K = find(m==0, 1 );
    else
        u1 = (u-p0);
        ind = find(cumsum(n/const)>=u1, 1 );        
        K = c(ind);
    end
end

function some_plot(z, U_mu, m, c, k, i, cmap)
        ind=find(m);
%         m(ind)
        hold off
%         length(ind)
        for j=1:length(ind)
            plot(z(1,c==ind(j)),z(2,c==ind(j)),'.','color',cmap(mod(5*ind(j),63)+1,:), 'markersize', 15);        
            hold on
            plot(U_mu(1,ind(j)),U_mu(2,ind(j)),'.','color',cmap(mod(5*ind(j),63)+1,:), 'markersize', 30); 
            plot(U_mu(1,ind(j)),U_mu(2,ind(j)),'ok', 'linewidth', 2, 'markersize', 10); 
        end
        plot(z(1,k),z(2,k),'or', 'linewidth', 3)            
        title(['i=' num2str(i) ',  k=' num2str(k) ', Nb of clusters: ' num2str(length(ind))]);
        xlabel('X')
        ylabel('Y') 
            xlim([-3 3]);
            ylim([-3 3]);
        pause(.01)
end
