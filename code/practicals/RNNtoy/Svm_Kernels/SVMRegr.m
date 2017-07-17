clc
clear all
close all
%Cp=textread('class_uno.txt');
Cp=textread('sinc.txt');
X=Cp(:,1:end-1);
Y=Cp(:,end);

    [N,D]  = size(X);
    Ntrain = round(N*0.7);
    Ntest  = N-Ntrain;
    index  = randperm(N);
    Xtrain = X(index(1:Ntrain),:);
    Xtest  = X(index(Ntrain+1:end),:);
    Ytrain = Y(index(1:Ntrain),:);
    Ytest  = Y(index(Ntrain+1:end),:);

[N,D]=size(Xtrain);
sigma2=0.7;
K=RBFcompute(Xtrain,sigma2);
C=10;
ep=0.05;
%% Variables para el Quadprog
K2=[K K
    K K];
un0=ones(N,1);
H=([un0;-un0]*[un0;-un0]').*K2;
eps=ep*ones(2*N,1);
Ym=[un0;-un0].*[Ytrain;Ytrain];
f=-(Ym-eps);
Aeq=[un0;-un0]';
Beq=0;
Lb=zeros(2*N,1);
Ub=C*ones(2*N,1);
options = optimset;
options=optimset(options,'Algorithm','interior-point-convex');
%options.MaxIter=1000;
An=quadprog(H,f,[],[],Aeq,Beq,Lb,Ub,[],options);

An(abs(An)<1e-4)=0;
[trash,sp]=find(An~=0);
Cvs=sum(sp);
an=An(1:N);
ahn=An(N+1:end);
[sv,sp]=find((an-ahn)~=0);
b=0;
Kpred=RBFcompute(Xtrain,sigma2,Xtest);
Kt=K(sv,sv);
yt=Ytrain(sv);
%suma=0;
b3=Ytrain-ep*ones(length(Ytrain),1)-K*(an-ahn);
b3=mean(b3);

for i = 1 : length(Ytrain)
   suma = 0;
    for j = 1 : length(Ytrain)
        suma = suma + (an(i)-ahn(i))*K(i,j);
    end
    b = b + Ytrain(i) - ep - suma;
end
b2=b/length(Ytrain);

Ypred=Kpred'*(an-ahn) + b2;
for n = 1 : length(Ytest)
    suma = 0;
    for i = 1:length(Ytrain)
    suma = suma + (an(i)-ahn(i))*Kpred(i,n);
    end
    Yt2(n) = suma + b2;
end

Y1=[Ytrain;Ytest];
Y2=[Ytrain;Ypred];
Y3=[Ytrain;Yt2'];
Xx=[Xtrain;Xtest];
[a1,b1]=sort(Xx);
Y1=Y1(b1);
Xx=Xx(b1);
Y2=Y2(b1);
[Sv1,Sp1] = find(an~=0);
[Sv2,Sp2] = find(ahn~=0);
Sv=sort([Sv1;Sv2]);
plot(Xx,Y2,'.k','markersize',12);
hold on
plot(Xx,Y1,'.r','markersize',12);
hold on
plot(Xx(Sv),Y1(Sv),'ok','markersize',18);
figure 
plot(abs((Ytest-Ypred)./Ytest))

