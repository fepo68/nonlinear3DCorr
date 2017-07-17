%% demRLS_regularización.m

clc
clear all
close all

x=linspace(0,10,200)';
y=3*sin(2*pi*x).*exp(-x/4);

%plot(x,y)
t=y+0.1*rand(200,1);

%% Entrenamiento 
N=40;
indexN=randperm(200);
indexN=sort(indexN(1:N)); %escoge 10 valores aleatorios entre 1 y 200
%y los organiza

xN=x(indexN); %datos de muestra para entrenamiento
tN=t(indexN);
%plot(xN,tN,'.b','markersize',20);

M=20;

muBasis=linspace(0,10,M-1)'; %se determinan las medias con separaciones iguales
sigmaBasis=0.8*(muBasis(2)-muBasis(1));%porcentaje de la diferencia entre dos mu continuos
phi=zeros (length(x),M-1);
for a=1:M-1
    %phi(:,a)=exp(-((x-muBasis(a)).^2)/(2*sigmaBasis^2));
    phi(:,a)=1./(1+exp(-((x-muBasis(a)))./(sigmaBasis)));
end
% figure
% plot(x,phi);

phiN=zeros(N,M);
phiN(:,1)=1;

for a=2:M
    %phiN(:,a)=exp(-((xN-muBasis(a-1)).^2)/(2*sigmaBasis^2));
    phiN(:,a)=1./(1+exp(-((xN-muBasis(a-1)))./(sigmaBasis)));
end
%se genera la matriz con las funciones base y con los datos de
%entrenamiento
lamda=0.1e-2;
Wml=(phiN'*phiN+lamda*eye(M))\phiN'*tN; %se encuentra la seudoinversa de la matriz phiN para encontrar el vector Wml

%% validación

phiT=[ones(length(x),1),phi]; %equivale a phi(x)
%matriz de tamaño: número de muestras por número de funciones base
yT=phiT*Wml;

figure;
plot(x,y,'r','linewidth',2);
hold on
plot(xN,tN,'.b','markersize',20);
hold on
plot(x,yT,'k','linewidth',2);