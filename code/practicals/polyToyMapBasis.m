function PHI_X = polyToyMapBasis(X,kern,D,fBasis)
PHI_X = {};
% fBasis = 'polynomial';
Mu1 = mean(X{1});
Mu2 = mean(X{2});
s = 0.01;


for d = 1:D
    Xd = X{d};
    [Nd,Md] = size(Xd)
    %     phiX = zeros(Nd(d),kern.degree);
    phiX =[];
    for n = 1:Nd
        x = Xd(n,:);
        
        
        
        if strcmp(fBasis,'polynomial')
%             xmap = designMatrix(x,'polynomial',kern,kern.degree);
            m = 0:kern.degree-1;
            xmap = kern.variance*(x*x'*kern.weightVariance+kern.biasVariance).^m;
        elseif strcmp(fBasis,'sigmoid')
            
            xmap = designMatrix(x,'sigmoid',Mu1,Mu2,s);
        end
        
        
        %         phiX(n,:) = xmap;
        phiX =[phiX;xmap];
    end
    
    PHI_X{d} = phiX;
end