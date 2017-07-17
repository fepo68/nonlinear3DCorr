function mySpyPlot(A,mrksize,Z1,Z2)
% function to plot an adjacency matrix given the clustering in Z
%
% A         I x J adjacency matrix
% Z         clustering assignment matrix
% mrksize   size of dots in spy-plot

if nargin<3
    Z1=[];
end
if nargin<4
    Z2=Z1;
end

if ~iscell(A)
   B=A;
   clear A;
   A{1}=B;
end
if nargin<2
    mrksize=1000/max(size(A{1}));
end
if isempty(mrksize)
    mrksize=1000/max(size(A{1}));
end
nn=length(A);
for n=1:length(A)
    if nn>1
       subplot(ceil(sqrt(nn)),ceil(sqrt(nn)),n);
       title(['A\{' num2str(n) '\}'],'FontWeight','bold')
    end
    [I,J]=find(A{n});
    linecol = [.8 .8 .8];
    hold on;
    [N1 N2]=size(A{n});
    if ~isempty(Z1)
        sZ1=cumsum(sum(Z1,2));
        sZ2=cumsum(sum(Z2,2));
        for k=1:length(sZ1)-1
           plot([0 sZ2(end)+1], N1+1-[sZ1(k)+.5 sZ1(k)+.5],'-','LineWidth',1,'Color',linecol); 
        end
        for k=1:length(sZ2)-1
           plot([sZ2(k)+.5 sZ2(k)+.5],N1+1-[0 sZ1(end)+1],'-','LineWidth',1,'Color',linecol);
        end
    end
    plot(J,N1+1-I,'.','MarkerSize',mrksize);
    plot([0 N2+1], [0 0],'-k','LineWidth',2);
    plot([0 0], [0 N1+1],'-k','LineWidth',2); 
    plot([N2+1 N2+1], [0 N1+1],'-k','LineWidth',2); 
    plot([0 N2+1], [N1+1 N1+1],'-k','LineWidth',2); 
    axis equal;
    axis tight;
    %h = xlabel(['N=' num2str(N) ', |{\bf A}|=' num2str(length(I)) ', L=' num2str(size(Z,1))]);
    set(gca,'XTick',[])
    set(gca,'XTickLabel',[])
    set(gca,'YTick',[])
    set(gca,'YTickLabel',[])
end