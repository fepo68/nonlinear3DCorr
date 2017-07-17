function ll = log_likelihoodEq4(params)

probNew = ((params.auxSumD)*log(2*pi))+(log(params.r)*(0.5*params.K*params.J))+((params.a*log(params.b))-...
                (params.ap*log(params.bp)))+gammaln(params.ap)-gammaln(params.a);
            
auxDetCj = 0;
for i = 1:params.J
    auxDetCj = auxDetCj+(log(det(inv(params.invCj(:,:,i))))*(0.5));
end

ll = probNew+auxDetCj;