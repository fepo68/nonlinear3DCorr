function Cov=RBFcompute(Xt,sigma,X2)
if nargin <3
    X2=Xt;
end
n2=dist2(Xt,X2);
Cov=exp(-0.5*sigma*n2); %% Se quita el 0.5 para el caso kernel!! 
