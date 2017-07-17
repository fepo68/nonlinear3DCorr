function Cov = Matern32compute(Xt,l,X2)
if nargin <3
    X2 = Xt;
end

if  l > 0
    n2 = dist2(Xt,X2);
    r = n2.^(1/2);
    Cov = (1 + sqrt(3)*n2/(2*l))*exp(-sqrt(3)*n2/(2*l));
    
else
    error(' L must be a Positive integer ');
end
