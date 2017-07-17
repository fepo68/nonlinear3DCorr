function Cov = RQcompute(Xt,l,alpha,X2)
if nargin <4
    X2 = Xt;
end
if alpha > 0 && l > 0
    
    n2 = dist2(Xt,X2);
    Cov = (1 + n2/(2*l*alpha)).^(-alpha);
else
    error(' Alpha and L must be Positive integers ');
end
