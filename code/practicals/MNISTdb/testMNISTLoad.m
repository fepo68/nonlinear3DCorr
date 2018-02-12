%% Test mnist
clear all
close all
clc

images = 'train-images-idx3-ubyte';
labels = 'train-labels-idx1-ubyte';
[ dataset ] = convertMNIST( images, labels );