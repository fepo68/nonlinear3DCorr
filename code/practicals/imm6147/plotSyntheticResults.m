function plotSyntheticResults(A,West,Ztrue,Zestimated)
% function to visualize the analysis of synthetically generated IRM data
%
% Usage:
%    plotSyntheticResults(A,Ztrue,Zestimated)
%
% Input:
%   A            Generated adjancency matrix
%   West         Link-predicted values
%   Ztrue        True assignment matrix
%   Zestimated   Estimated assignment matrix
%
% Written by Morten Mørup

J=size(A,1);
figure;
subplot(2,3,1);
mySpyPlot(A,1000/J);
title('Generated Graph','FontWeight','Bold')
subplot(2,3,2);
[val,ind]=sort(sum(Ztrue,2),'descend');
Ztrue=Ztrue(ind,:);
[A_sorted,Z_True_sorted]=sortGraphUnipartite(A,Ztrue,rand(size(Ztrue,1)));
mySpyPlot(A_sorted,1000/J,Z_True_sorted);
title('Correctly Sorted Generated Graph','FontWeight','Bold')
subplot(2,3,3);
[A_sorted,Z_Est_sorted]=sortGraphUnipartite(A,Zestimated,rand(size(Zestimated,1)));
mySpyPlot(A_sorted,1000/J,Z_Est_sorted);
title('Estimated Sorting of the Generated Graph','FontWeight','Bold')
subplot(2,3,4);
imagesc(Z_True_sorted*Z_Est_sorted'); colormap(1-gray); axis equal; axis tight; title(['Ztrue*Zestimated^T, NMI=' num2str(round(calcNMI(Z_True_sorted,Z_Est_sorted)*100)/100) ', AUC=' num2str(round(calcAUC(West,A)*100)/100)],'FontWeight','Bold')
subplot(2,3,5);
imagesc(Z_True_sorted); colormap(1-gray); axis off; title('True (sorted) Assigment Matrix Z','FontWeight','Bold')
subplot(2,3,6);
imagesc(Z_Est_sorted); colormap(1-gray); axis off; title('Estimated (sorted) Assignment Matrix Z','FontWeight','Bold')
