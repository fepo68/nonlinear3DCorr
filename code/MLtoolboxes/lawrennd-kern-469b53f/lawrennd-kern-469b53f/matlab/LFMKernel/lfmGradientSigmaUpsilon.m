function    g = lfmGradientSigmaUpsilon(gamma, sigma2, Tt1, Tt2);

% LFMGRADIENTSIGMAUPSILON Gradient of the function \upsilon(z) with respect
% to \sigma.
% FORMAT
% DESC Computes the gradient of the function \upsilon(z) with respect to
% the length-scale of the input "force", \sigma.
% ARG gamma : Gamma value of the system.
% ARG sigma2 : length scale of latent process.
% ARG Tt1 : first time input (number of time points 1 x number of time points 2).
% ARG Tt2 : second time input (number of time points 1 x number of time points 2).
% RETURN g : Gradient of the kernel with respect to the desired
% parameter.
%
% COPYRIGHT : David Luengo, 2008
%
% SEEALSO : lfmKernGradient, lfmXlfmKernGradient, lfmGradientSigmaH

% LFM


%Parameters of the function

sigma = sqrt(sigma2);

% Initialization of vectors and matrices

g = zeros(size(Tt1));

Z1 = (Tt1-Tt2)/sigma - sigma*gamma/2;
Z2 = Tt2/sigma + sigma*gamma/2;

%%% Gradient Evaluation %%%

% Evaluation of the gradient when real(Z1)>=0 and real(Z2)>=0

ind = (real(Z1)>=0) & (real(Z2)>=0);
if any(any(ind))
    g(ind) = sigma*(gamma^2)*exp(sigma2*(gamma^2)/4 ...
            - gamma*(Tt1(ind)-Tt2(ind))) ...
        - 2*exp(-(Tt1(ind)-Tt2(ind)).*(Tt1(ind)-Tt2(ind))/sigma2 ...
            + log((Tt1(ind)-Tt2(ind)).*(Tt1(ind)-Tt2(ind)).*W(j*Z1(ind))/(sigma^3) ...
            + ((Tt1(ind)-Tt2(ind))/sigma2+gamma/2) ...
            .* (1/sqrt(pi)-Z1(ind).*W(j*Z1(ind))))) ...
        - 2*exp(-Tt2(ind).*Tt2(ind)/sigma2-gamma*Tt1(ind) ...
            + log(Tt2(ind).*Tt2(ind).*W(j*Z2(ind))/(sigma^3) ...
            + (Tt2(ind)/sigma2-gamma/2) .* (1/sqrt(pi)-Z2(ind).*W(j*Z2(ind)))));
end

% Evaluation of the gradient when real(Z1)<0 and real(Z2)>=0

ind = (real(Z1)<0) & (real(Z2)>=0);
if any(any(ind))
    g(ind) = 2*exp(-(Tt1(ind)-Tt2(ind)).*(Tt1(ind)-Tt2(ind))/sigma2 ...
            + log((Tt1(ind)-Tt2(ind)).*(Tt1(ind)-Tt2(ind)).*W(-j*Z1(ind))/(sigma^3) ...
            - ((Tt1(ind)-Tt2(ind))/sigma2+gamma/2) ...
            .* (1/sqrt(pi)+Z1(ind).*W(-j*Z1(ind))))) ...
        - 2*exp(-Tt2(ind).*Tt2(ind)/sigma2-gamma*Tt1(ind) ...
            + log(Tt2(ind).*Tt2(ind).*W(j*Z2(ind))/(sigma^3) ...
            + (Tt2(ind)/sigma2-gamma/2) .* (1/sqrt(pi)-Z2(ind).*W(j*Z2(ind)))));
end

% Evaluation of the rest of the gradient when real(Z1)>=0 and real(Z2)<0

ind = (real(Z1)>=0) & (real(Z2)<0);
if any(any(ind))
    g(ind) = - 2*exp(-(Tt1(ind)-Tt2(ind)).*(Tt1(ind)-Tt2(ind))/sigma2 ...
            + log((Tt1(ind)-Tt2(ind)).*(Tt1(ind)-Tt2(ind)).*W(j*Z1(ind))/(sigma^3) ...
            + ((Tt1(ind)-Tt2(ind))/sigma2+gamma/2) ...
            .* (1/sqrt(pi)-Z1(ind).*W(j*Z1(ind))))) ...
        + 2*exp(-Tt2(ind).*Tt2(ind)/sigma2-gamma*Tt1(ind) ...
            + log(Tt2(ind).*Tt2(ind).*W(-j*Z2(ind))/(sigma^3) ...
            - (Tt2(ind)/sigma2-gamma/2) .* (1/sqrt(pi)+Z2(ind).*W(-j*Z2(ind)))));
end

% Evaluation of the rest of the gradient when real(Z1)>=0 and real(Z2)>=0

ind = (real(Z1)<0) & (real(Z2)<0);
if any(any(ind))
    g(ind) = - sigma*(gamma^2)*exp(sigma2*(gamma^2)/4 ...
            - gamma*(Tt1(ind)-Tt2(ind))) ...
        + 2*exp(-(Tt1(ind)-Tt2(ind)).*(Tt1(ind)-Tt2(ind))/sigma2 ...
            + log((Tt1(ind)-Tt2(ind)).*(Tt1(ind)-Tt2(ind)).*W(-j*Z1(ind))/(sigma^3) ...
            - ((Tt1(ind)-Tt2(ind))/sigma2+gamma/2) ...
            .* (1/sqrt(pi)+Z1(ind).*W(-j*Z1(ind))))) ...
        + 2*exp(-Tt2(ind).*Tt2(ind)/sigma2-gamma*Tt1(ind) ...
            + log(Tt2(ind).*Tt2(ind).*W(-j*Z2(ind))/(sigma^3) ...
            - (Tt2(ind)/sigma2-gamma/2) .* (1/sqrt(pi)+Z2(ind).*W(-j*Z2(ind)))));
end
