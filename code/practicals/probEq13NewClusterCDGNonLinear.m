function [probNew,ap_sdn_j,bp_sdn_j, mu_j_sdn_j, invCj_sdn_j] = probEq13NewClusterCDGNonLinear(params,W,X,d,n)

% [apNew,bpNew,mu_jNew,invCjNew] = equation4IWataNewCluster_no_dn(params,W,X,d,n);

xdn = X{d}(n,:)';

ap_sdn_j = params.ap; % eq 14
% for C_j
kern = kernCreate(W{d}'*W{d},params.kernType);
kern.variance = params.variance;
kern.inversewidth = params.inversewith;
Kd = kernCompute(kern,W{d}');
% invCj_sdn_j = W{d}'*W{d} + params.r*eye(params.K); % eq 17
invCj_sdn_j = Kd + params.r*eye(params.K); % eq 17
% for mu
kern = kernCreate(W{d}',params.kernType);
kern.variance = params.variance;
kern.inversewidth = params.inversewith;
Kwx = kernCompute(kern,W{d}',xdn');
% mu_j_sdn_j = invCj_sdn_j\(W{d}'*xdn);
mu_j_sdn_j = invCj_sdn_j\(Kwx);
% for bp
kern = kernCreate(xdn',params.kernType);
kern.variance = params.variance;
kern.inversewidth = params.inversewith;
kxx = kernCompute(kern,xdn');
% bp_sdn_j = params.bp_no_dn +0.5*(xdn'*xdn)-0.5*mu_j_sdn_j'*invCj_sdn_j*mu_j_sdn_j;
bp_sdn_j = params.bp_no_dn +0.5*(kxx)-0.5*mu_j_sdn_j'*invCj_sdn_j*mu_j_sdn_j;


% probNew = ((2*pi)^(-0.5*params.Md(d)))*(params.r^(0.5))*exp(((params.ap_no_dn*log(params.bp_no_dn))-...
%     (ap_sdn_j*log(bp_sdn_j))))*exp(gammaln(ap_sdn_j)-gammaln(params.ap_no_dn))*...
%     ((1/det(invCj_sdn_j))^0.5); % without rI


probNew = ((2*pi)^(-0.5*params.Md(d)))*(params.r^(0.5))*exp(((params.ap_no_dn*log(params.bp_no_dn))-...
    (ap_sdn_j*log(bp_sdn_j))))*exp(gammaln(ap_sdn_j)-gammaln(params.ap_no_dn))*...
    ((1/det(invCj_sdn_j))^0.5)/((1/det(params.r*eye(params.K)))^0.5);