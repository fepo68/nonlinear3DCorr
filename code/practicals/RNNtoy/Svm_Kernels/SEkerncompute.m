function Cov=SEkerncompute(Xt,sigma2,l,X2)
if nargin <4
    X2=Xt;
end

n2=dist2(Xt,X2);
Cov=sigma2*exp(-0.5*n2/(l*l)); 
