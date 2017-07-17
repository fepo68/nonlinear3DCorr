%% Databases to UCM-LVM model
clear all
close all
clc

% load  dataToyIwatamodel27Ene162D40Ndj5J10K

dataSets = {'iris','wine','glass'};
dataSet = dataSets{2};
if dataSet == 'iris'
    load iris_dataset
    D = 2;
    Nd = [150,150];
    Md = [2,2];
    X{1} = irisInputs(1:2,:)';
    X{2} = irisInputs(3:4,:)';
    J = 3;
    K = 3;
    W{1} = eye(Md(1),K);
    W{2} = eye(Md(2),K);
    
    S{1} = [ones(1,50),2*ones(1,50),3*ones(1,50)];
    S{2} = [ones(1,50),2*ones(1,50),3*ones(1,50)];
end
if dataSet =='wine'
    load wine_dataset
    [targets,~] = find(wineTargets==1);
    S{1} = targets';
    S{2} = targets';
    id1 = randi(13,[1,6]);
    X{1} = wineInputs(id1,:)';
    X{2} = wineInputs(setdiff(1:13,id1),:)';
    D = 2;
    Nd = [178,178];
    Md = [6,7];
    J = max(targets);
    K = 3;
    W{1} = eye(Md(1),K);
    W{2} = eye(Md(2),K);
end
save(['dataToy',dataSet,'Dataset']);