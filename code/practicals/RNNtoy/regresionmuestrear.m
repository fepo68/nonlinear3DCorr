clc
clear all
close all

load ejemplo_class_dos.mat
N=30;
index=sort(randperm(size(X,1),N));
XTrain=X(index(1:N),1); %tomar los datos de entrenamiento
XTrain=[ones(N,1) XTrain]; %sumarles la columna de 1

tTrain=t(index(1:N)); %tomar los t de entrenamiento correspondientes a los X tomados (el mismo index)
%% Escoger datos TEST
XTest=[ones(length(t)-N,1) X(index(N+1:end),:)];
tTest=t(index(N+1:end));

W=zeros(size(XTrain,2),1); %se inicializa el W (dimensionalidad dada x las filas)


%% Muestra del Entrenamiento
gradv(gradv==-1)=[];%eliminar los valores que son negativos
numiters=1:length(gradv);
figure
plot(numiters,gradv,'k','linewidth',1.5);
xlabel('Numero de Iteraciones');
ylabel('Norma Cuadrada del gradiente');
grid on
hold on
figure
plot(XTrain(tTrain==0,2),XTrain(tTrain==0,3),'.b','markersize',10);
hold on
plot(XTrain(tTrain==1,2),XTrain(tTrain==1,3),'xr','linewidth',2);
grid on



%% Trazado Recta de Carga

X1R=linspace(-3,3,100)';
X2R=-(W(2)./W(3))*X1R-(W(1)./W(3));
%% GRAFICA
t2=(t+1)/2; %transforma los -1 en ceros
N=10;
%% ESCOGER DATOS Train
index=randperm(length(t)); %indexar valores


plot(X1R,X2R,'k','linewidth',3.2);
grid on

%% Validacion con los Test

 a=XTest*W;
 yTest=1./(1+exp(-a));
 yTestP=zeros(length(yTest),1); %ya estan los ceros 
 yTestP(yTest>0.5)=1; %se ponen en 1 las que son mayores a 0.5
 
 error=100*sum(abs(yTestP-tTest))/length(tTest);
