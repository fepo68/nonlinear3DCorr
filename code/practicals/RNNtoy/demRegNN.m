%DEMREGRNN

clc
clear all

close all
D=1;
M=10;
K=1;

sizew1=(D+1)*(M+1);
sizew2=(M+1)*K;
sizew=sizew1+sizew2;


%% generacion del vector de pesos

alpha=0.1;
w=gsamp(randn(sizew,1),alpha*eye(sizew),11)'; % gausiana multivariada, con media rand sizew y desviacion standard con alpha 0.o

w1=w(1:sizew1); %se genera vector
W1=reshape(w1,M+1,D+1);
w2=w(sizew1+1:sizew);
W2=reshape(w2,K,M+1);

%% Definicion del espacio de entrada
N=100;
x=linspace(-5,5,N)'; %vector de entrada X (con 100 elementos x defecto)
% el X que se necesita será  con el bias de 1 nos

X=[ones(N,1) x];

%existen j combinaciones lineales (M ) de los datos de entrada 

aj=W1*X'; %activaciones de tamanyo M+1 x N por tamanyo

zj=tanh(aj); %combinacion no lineal de las entradas ,funciones de las activaciones! %M+1 x N % funciones de las activaciones o funciones base

ak=W2*zj; %activaciones de la segunda capa tamanyo K*N

% como es un modelo de regresion yk = ak

yk=ak'; %obviamente es una combinacion lineal de todas esas tangentes

%% propagacion hacia adelante finished, now plot

plot(x,yk,'k','linewidth',1.5);
hold on
plot(x,zj') ;

