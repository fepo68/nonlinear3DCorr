%DEM Machine Learning Reg1
clc
clear
close all
M=5; %funcion base
N=100; %numero de puntos
x=linspace (-1,1,N)'; %comando de entrada (la comilla es para transponer ya que matlab todo me lo muestra columna)

%polinomial

 fi= zeros(N,M-1);
for i=1:M-1
   
    fi(:,i)=x.^i;%cuando es elemento a 
                 %elemento pongo el punto, aqui almaceno los datos
                 %preguntar lo de los 2 puntos
end
subplot(3,1,1)% columna,fila,posicion
plot(x,fi)
title('Polinomial')
grid

%exponencial

fi(:,i)=x.^i;
miu=linspace(-1,1,M-1);
s2=0.1;% desviacion standar
 for i=1:M-1
fi(:,i)=exp(-(x-miu(i)).^2/(2*s2));
 end

 subplot(3,1,2)
 plot(x,fi)
 title('Exponencial')
 
 % Sigmoidal
 
fi = zeros (N,M-1);

miu = linspace(-1,1,M-1);
s=0.1;
for i=1:M-1
fi(:,i)=1./(1+exp(-(x-miu(i))/s));
end
subplot(3,1,3)
plot(x,fi)
title('Sigmoidal')
grid

