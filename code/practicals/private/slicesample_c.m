function [c, m, U_mu, U_Sigma, U_SS] = slicesample_c(y, c, U_mu, U_Sigma, U_SS, hyperG0, alpha)

% Sample the allocation variables
% Use slice sampling, cf. (Kalli, Griffin & Walker)

[~, n] = size(y);
MAX_CLUST = size(U_mu, 2);

% MAX_CLUST = n;
ind = unique(c);
N = length(ind);
m = zeros(MAX_CLUST, 1);
for t=1:N
    m(ind(t)) = sum(c==ind(t));    
end

% Sample the weights
pi = zeros(MAX_CLUST, 1);
temp = gamrnd([m(ind); alpha], 1);
temp = temp/sum(temp);
pi(ind) = temp(1:end-1);
R  = temp(end);

% Sample the latent u
u = rand(n,1).*pi(c);

% Sample the remaining weights that are needed with stick-breaking
ind_new = find(m==0);
t=0;
while((sum(pi)<1-min(u)) && (t<length(ind_new)))
    t = t+1;
    beta_temp = betarnd(1, alpha);
    pi(ind_new(t)) = R*beta_temp;
    R = R * (1-beta_temp);
    % Sample new atoms
    [U_mu(:, ind_new(t)), U_Sigma(:, :, ind_new(t))] = normalinvwishrnd(hyperG0);
end
% t
% pause
ind_new = ind_new(1:t);


% Get the likelihood 
ind_all = [ind; ind_new];
for i=1:n
%     likelihood(i,:) = mvnpdf(repmat(y(:,i), 1, length(ind_all))', U_mu(:, ind_all)', U_Sigma(:, :, ind_all))';
    loglikelihood(i,:) = loggausspdf(repmat(y(:,i), 1, length(ind_all))', U_mu(:, ind_all)', U_Sigma(:, :, ind_all))';
end
loglikelihood = loglikelihood - max(loglikelihood(:));
likelihood = exp(loglikelihood);
% likelihood
% pause

[A, B] = meshgrid(pi(ind_all), u);
proba = (A>B) .* likelihood;
for i=1:n
    ind1 = find(rand*sum(proba(i, :))<=cumsum(proba(i, :)), 1);
    if isempty(ind1)
        keyboard
    end
    m(c(i)) = m(c(i)) - 1;
    U_SS(c(i)) = downdate_SS(y(:,i),U_SS(c(i)));
    
    c(i) = ind_all(ind1);
    U_SS(c(i)) = update_SS(y(:,i), U_SS(c(i)));
    m(c(i)) = m(c(i)) + 1;
end