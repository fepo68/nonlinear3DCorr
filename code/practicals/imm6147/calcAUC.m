function [AUC,TPR,FPR]=calcAUC(West,A)
% Calculates the Area Under Curve (AUC) of the Receiver Operating
% Characteristic (ROC)
%
% Usage:
%   [AUC,TPR,FPR]=calcAUC(West,A)
% 
% Input: 
%   West:   Graph(s) of estimated probabilities of generating links in the graph for
%           entries treated as missing
%   A:      True Graph(s)
%
% Output:
%   AUC: The area under the ROC curve
%   TPR: True Positive Rate
%   FPR: False Positive Rate
%
% Written by Morten Mørup

if iscell(A)
    MAP=[];
    C=[];
    for n=1:length(A)
        [a,b,map]=find(West{n});
        MAP=[MAP; map];
        Wlogical=West{n}>0;
        [a,b,c]=find(Wlogical.*A{n}+Wlogical);
        C=[C; c-1];
    end
else
    [a,b,MAP]=find(West{1});
     Wlogical=West{1}>0;
    [a,b,C]=find(Wlogical.*A+Wlogical);
    C=C-1;
end
if ~isempty(MAP)
    [val,ind]=sort(MAP(:),'ascend');
    x=C(ind);
    N0=sum(1-x);
    N1=sum(x);
    FNR=[zeros(length(x),1); 1];
    TNR=[zeros(length(x),1); 1];
    for k=1:length(val)
        ind=find(val<val(k));    
        FNR(k+1)=sum(x(ind))/N1;
        TNR(k+1)=sum(1-x(ind))/N0;
    end

    TPR = 1-FNR;
    FPR = 1-TNR;
    AUC = -trapz(FPR,TPR);
else
    AUC=NaN;
    TPR=NaN;
    FPR=NaN;
end


