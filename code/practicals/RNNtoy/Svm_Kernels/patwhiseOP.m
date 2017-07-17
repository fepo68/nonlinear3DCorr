clc 
clear all
close all


Cp=textread('diabetes.txt');

seed=2240; 
% seed=1320; 
rand('seed',seed);
randn('seed',seed);

%% Normalizacion
X=Cp(:,1:end-1);
Y=Cp(:,end);
[N,D]=size(X);

%normalizacion X
varX=zeros(1,D);
mX=sum(X,1)/N;
for i=1:N
   varX=varX+(X(i,:)-mX).^2;
end
varX=varX/N;
desX=sqrt(varX);
X=(X-repmat(mX,N,1))./(repmat(desX,N,1));

%normalizacion Y
mY=sum(Y)/length(Y);
Y=Y-mY;

%check in  OK 

%% Training

% Se eligen la cantidad de elementos de entrenamiento y validacion
Ntrain=round(N*0.7);
Ntest=N-Ntrain;
index=randperm(N);
Xtrain=X(index(1:Ntrain),:);
Ytrain=Y(index(1:Ntrain),:);
Xtest=X(index(Ntrain+1:end),:);
Ytest=Y(index(Ntrain+1:end),:);


%% lambda maximo y minimo
eps=0.001;
K=100;
Xk=Xtrain'*Ytrain;
ldamax=max(Xk/Ntrain);
ldamin=eps*ldamax;

% hallar los lambdas vectores y graficar el regularization path
W=CoordinateDescent(Xtrain,Ytrain,ldamax);

Lk=linspace(ldamax,ldamin,K);
Wk=zeros(K,size(X,2));
for i=1:K
  Wz=CoordinateDescent(Xtrain,Ytrain,Lk(i));
  Wk(i,:)=Wz;
  Cte=Lk(i)*norm(Wk(i,:),1);
  Rss=sum((Ytest-Xtest*Wk(i,:)').^2)/(2*Ntest);
  E(:,i)=Rss+Cte;
%   if i>3
%       break
%   end
end
figure (1)
plot(Wk(1:end,:))
%% Ridge regression 
Wk2=zeros(K,10);
for i=1:K
  Wk2(i,:)=(Xtrain'*Xtrain+Lk(i)*eye(D))^-1*(Xtrain'*Ytrain);  
  Cte=Lk(i)*norm(Wk2(i,:),1);
  Rss=sum((Ytest-Xtest*Wk2(i,:)').^2)/(2*Ntest);
  E2(:,i)=Rss+Cte;
end
figure (2)
plot(Wk2(1:end,:))





