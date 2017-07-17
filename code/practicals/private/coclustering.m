function [out, subindex] = coclustering(c, pas, c_hat)

%
% COCLUSTERING computes the coclustering, or similarity matrix
%   [out, subindex] = COCLUSTERING(c, pas, c_hat)
%
%   INPUTS
%   - c:    Matrix of size L*n. c(i,j) is the cluster allocation of object
%           i at iteration j
%   - pas:  Provides the coclustering every pas objects (Default=1)
%   - c_hat:Partition used to reorder the objects. 
%           Objects are reordered given the size of the cluster they belong to (bigger clusters first) 
%          (Default:c_hat=c(:,end))
%
%   OUTPUTS:
%   - out:  coclustering matrix. out(i,j) is the number of iterations where
%           c(i)=c(j)
%   - subindex: indices of the objects in the coclustering matrix
%
%   See also CLUSTER_EST_BINDER
%--------------------------------------------------------------------------
% EXAMPLE
% c_st=[1,1,2;1,2,3;2,2,1]
% coclust = coclustering(c_st);
%

% Copyright Francois Caron, University of Oxford, 2014
%--------------------------------------------------------------------------


if nargin<2
    pas = 1;
end
if nargin<3
    c_hat = c(:,end); % sorts the coclustering matrix by clusters of decreasing size
end

% Compute matrix of co-clustering from draws from the posterior partition

[n, ~] = size(c);
ind = (unique(c_hat));
N = length(ind);
 

%  m = zeros(50, size(lambda_st,3));
for t=1:N
    m(ind(t)) = sum(c_hat==ind(t)); 
end
[~, ind2] = sort(m, 'descend');
ind3 = zeros(n, 1);
k = 1;
for t=1:length(ind2)
    indt = find(c_hat==ind2(t));
    for i=1:length(indt)
        ind3(k) = indt(i);
        k=k+1;
    end
end


nmax = n;
subindex = 1:pas:nmax;
c2 = c(ind3(subindex), :);

n2 = length(subindex);
out = zeros(n2,n2,'int16');
for i=1:n2
    for j=i+1:n2        
        out(i, j) = sum(c2(i,:) == c2(j,:));
%         out(j, i) = out(i,j);
    end
end
out = out + out' + size(c, 2)*eye(size(out), 'int16');

