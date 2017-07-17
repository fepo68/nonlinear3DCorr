function [x,fval,exitflag,output,lambda] = AAA
%% This is an auto generated MATLAB file from Optimization Tool.

%% Start with the default options
options = optimset;
%% Modify options setting
options = optimset(options,'Display', 'off');
options = optimset(options,'Algorithm', 'interior-point-convex');
[x,fval,exitflag,output,lambda] = ...
quadprog([],[],[],[],[],[],[],[],[],options);
