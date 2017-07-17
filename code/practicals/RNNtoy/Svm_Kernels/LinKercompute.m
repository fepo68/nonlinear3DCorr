function Cov = LinKercompute(Xt,sigma1,sigma2,c,X2)
if nargin <5
    X2 = Xt;
end

dim1 = size(Xt,1);
dim2 = size(X2,1);

Cov = zeros (dim1,dim2);
if size(Xt,2) == size(X2,2)
    for i = 1 :dim1
        
        for j = 1 : dim2
            Cov(i,j) = sigma1^2 + (Xt(i,:) - c)*(X2(j,:) - c)'*sigma2^2 ;
        end
        
    end
else
    error('The Columns must be the same ');
end

