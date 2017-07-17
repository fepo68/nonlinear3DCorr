function logpdf = loggausspdf(x, mu, Sigma)

[n, d] = size(x);

[n2,d2] = size(mu);
% Check dimension and center data
if d2 ~= d 
    error(message('X and mu must have the same nb of columns'));
elseif n2 == n 
    x0 = x - mu;
elseif n2 == 1 % mean is a single row, rep it out to match data
    x0 = bsxfun(@minus,x,mu);
else % sizes don't match
    error(message('mu must has 1 row or the same nb as x'));
end



if ndims(Sigma)==2
    [R,err] = chol(Sigma);
    % Standardized data
    xstand =  x0 / R;
    if err~= 0
        error('Sigma must be semi-definite positive');
    end
    % log(sqrt(Sigma))
    logSqrtDetSigma = sum(log(diag(R)));
    
elseif ndims(Sigma)==3
    xstand = zeros(n,d);
    logSqrtDetSigma = zeros(n,1);
    for i = 1:n
        [R,err] = chol(Sigma(:,:,i));
        if err ~= 0
            error(message('stats:mvnpdf:BadMatrixSigmaMultiple'));
        end
        xstand(i,:) = x0(i,:) / R;
        logSqrtDetSigma(i) = sum(log(diag(R)));
    end
end

% Quadratic form
xquad = sum(xstand.^2, 2);
logpdf = -0.5*xquad - logSqrtDetSigma - d*log(2*pi)/2;