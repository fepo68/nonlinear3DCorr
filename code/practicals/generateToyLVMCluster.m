%% Genera Toy PhD. UCM-LVM, as in Iwata's Paper from 2016

clear all
close all
clc

J = 5;
D = 2;
K = 10;
Md = [50 50];
Nd = [200 200];

for D = 2:10 % for several domains experiments
    
    Md = 50*ones(1,D);
    Nd = 200*ones(1,D);
    
    for i = 1:5
        
        
        s = RandStream('mt19937ar','Seed',10e5*i);
        RandStream.setGlobalStream(s);
        
        
        [S,W,X] = createToyLVMCluster(Nd,Md,D,J,K);
        
        plot(X{1}(:,1),X{1}(:,2),'*r')
        
        save(['synth',num2str(K),'exp',num2str(i),'D',num2str(D),'.mat'],'S','W','X','Nd','Md','D','J','K');
    end
end