clear all
close all
clc
addpath(genpath('./starter'));
% Change the filenames if you've saved the files under different names
% On some platforms, the files might be saved as
% train-images.idx3-ubyte / train-labels.idx1-ubyte
images = loadMNISTImages('train-images-idx3-ubyte');
labels = loadMNISTLabels('train-labels-idx1-ubyte');

% We are using display_network from the autoencoder code
% display_network(images(:,100:200)); % Show the first 100 images
disp(labels(1:10));

%% Processing MNITS
% First we downsampled each image to have 16x16 size
for k = 1:400:2000
    X = {};
    S = {};
    Nd = [200,200];
    
    dataDs = images(:,k:k+399);
    [Md,N] = size(dataDs);
    dataDlowres = zeros(16^2,N);
    for i = 1:N
        img = reshape(dataDs(:,i),[28,28]);
        img = imresize(img,[16,16]);
        %     imshow(img,[]);
        dataDlowres(:,i) = img(:);
    end
    labelsDb = labels(1:2*Nd(1))+1;
    Md = [16^2,16^2];
    X{1} = dataDlowres(:,1:200)';
    X{2} = dataDlowres(:,201:end)';
    S{1} = labelsDb(1:200)';
    S{2} = labelsDb(201:end)';
    D = 2;
    save(['dataToy','MNIST','Dataset',num2str(k),'exp'],'D','X','S','Md','Nd');
    
end



