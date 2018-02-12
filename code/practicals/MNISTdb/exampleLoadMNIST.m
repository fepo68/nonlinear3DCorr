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
% First we downsampled each image to have pixNxpixN size
iRes = false;
count = 1;
for k = 1:800:4000
    X = {};
    S = {};
    Nd = [400,400];
    dataDs = 255*images(:,k:k+799);
    labelsDb = labels(k:k+799);
    [aa pp] = sort(labelsDb);
    dataDs = dataDs(:,pp);
    labelsDb = aa+1;
    CVO = cvpartition(labelsDb,'k',2);
    d1Idx = CVO.training(1);
    d2Idx = CVO.test(1);
    
    
    
    
    [Md,N] = size(dataDs);
    pixN = 28;
    dataDlowres = zeros(pixN^2,N);
    if iRes ==true
        for i = 1:N
            img = reshape(dataDs(:,i),[28,28]);
            %         img = imresize(img,[pixN,pixN]);
            %     imshow(img,[]);
            dataDlowres(:,i) = img(:);
        end
    else
        dataDlowres = dataDs;
    end
    %     labelsDb = labels(1:2*Nd(1))+1;
    
    Md = [pixN^2,pixN^2];
    X{1} = dataDlowres(:,d1Idx)';
    X{2} = dataDlowres(:,d2Idx)';
    S{1} = labelsDb(d1Idx)';
    S{2} = labelsDb(d2Idx)';
    D = 2;
    save(['dataToy','MNIST','Dataset',num2str(count),'exp'],'D','X','S','Md','Nd');
    count = count +1;
    
end



