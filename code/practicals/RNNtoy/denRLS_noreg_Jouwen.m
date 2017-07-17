%% demRLS.m

clc
clear all
close all

x=linspace(0,1,200)';
y=sin(2*pi*x);

%plot(x,y)
t=y+0.1*rand(200,1);

%% Entrenamiento 
N=30;
indexN=randperm(200);
indexN=sort(indexN(1:N)); %escoge 10 valores aleatorios entre 1 y 200
%y los organiza

xN=x(indexN); %datos de muestra para entrenamiento
tN=t(indexN);
%plot(xN,tN,'.b','markersize',20);

M=10;

muBasis=linspace(0,1,M-1)'; %se determinan las medias con separaciones iguales
sigmaBasis=0.5*(muBasis(2)-muBasis(1));%porcentaje de la diferencia entre dos mu continuos
phi=zeros (length(x),M-1);
for a=1:M-1
    phi(:,a)=exp(-((x-muBasis(a)).^2)/(2*sigmaBasis^2));
end

%plot(x,phi);

phiN=zeros(N,M);
phiN(:,1)=1;

for a=2:M
    phiN(:,a)=exp(-((xN-muBasis(a-1)).^2)/(2*sigmaBasis^2));
end
%se genera la matriz con las funciones base y con los datos de
%entrenamiento
Wml=(phiN'*phiN)\phiN'*tN; %se encuentra la seudoinversa de la matriz phiN para encontrar el vector Wml

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