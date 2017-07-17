function kern = lfmKernParamInit(kern)

% LFMKERNPARAMINIT LFM kernel parameter initialisation.
% The latent force model (LFM) kernel is the result of a second order
% differential equation where there is assumed to be a force driving the
% system which is drawn from a Gaussian process with an RBF kernel,
% i.e. we have the following differential equation,
%
% S f(t-delta) = m x''(t) + cx'(t) + kx(t),
%
% where m is a mass, c is a damping coefficient, k is a spring constant,
% S is a scalar sensitivity and delta is a time delay
% and B is an initial level. 
  
% If f(t) is assumed to come from a
% Gaussian process with an RBF covariance function x(t) is a Gaussian
% process with a covariance function provided by the single latent force
% model kernel.
%
% The kernel is designed to interoperate with the multiple output
% block kernel so that f(t) can be inferred given several different
% instantiations of x(t).
%
% The parameters (m, c, delta and k) are constrained positive.
%
% FORMAT
% DESC initialises the latent force model kernel structure with some
% default parameters.
% ARG kern : the kernel structure which requires initialisation.
% RETURN kern : the kernel structure with the default parameters placed in.
%
% SEEALSO : kernCreate, kernParamInit, lfmKernCompute
%
% COPYRIGHT : Neil D. Lawrence, 2007

% LFM

if kern.inputDimension > 1
  error('LFM kernel only valid for one-D input.')
end

kern.delay = 0;
kern.mass = 1;
kern.spring = 1;
kern.damper = 1;
kern.sensitivity = 1;

kern.initVal = 1;

kern.variance = 1;
kern.inverseWidth = 1;
kern.nParams = 5;

kern.transforms.index = [1 2 3 4];
kern.transforms.type = optimiDefaultConstraint('positive');

kern.isStationary = false;

% Force any precomputation contained in lfmKernExpandParam
params = lfmKernExtractParam(kern);
kern = lfmKernExpandParam(kern, params);
