function [Mclust,Xind,Yind,p,lx,ly]=myCluster(M,minV,diag,clusterRows,clusterColumns,varargin)

% Mclust=myCluster(M,minV,diag,cluster_rows,cluster_cols,corr_type,make_heatmap)
%
% e.g. myCluster(M,10,0,'euclidean','no')
% diag = 0/1
% corr_type: any of the offered types in standard pdist
% M needs to contain: M.data, M.rowlabels, M.collabels
% NaNs are converted to zeros before clustering
% will cluster and diagonalize data in M.data
% use minV for dumping data, will dump all rows and columns with max < minV

if isempty(varargin)
    corr_type='euclidean';
else
    corr_type=varargin{1};
end

%ns=sum(isnan(M.data),2);
%mn=b-min_non_nan;
%x=ns<=mn;
%fprintf([num2str(sum(x)) ' rows met your cuttoff of ' num2str(min_non_nan) ' non-NaNs' '\n']);
%M.data=M.data(x,:);
%M.rowlabels=M.rowlabels(x);

% subset data

y=max(abs(M.data))>minV;
x=max(abs(M.data'))>minV;

fprintf([num2str(sum(x)) ' rows met your cuttoff of ' num2str(minV) '\n']);
fprintf([num2str(sum(y)) ' columns met your cuttoff of ' num2str(minV) '\n']);

M.data=M.data(x,y);
M.rowlabels=M.rowlabels(x);
M.collabels=M.collabels(y);

% store NaNs but replace with zeros for clustering
nn=isnan(M.data);
M.data(nn)=0;

[a,b]=size(M.data);

% get Xind
lx=0;
p=0;
if clusterRows==1
    p=pdist(M.data,corr_type);
    lx=mylinkage(p,3);
    %[~,Xind,xx]=dendrogram(l,a);
    [~,Xind]=mydendrogram(lx,a);
    print -depsc dend.x.eps
else
    Xind=1:a;
end
    
% get Yind
ly=0;
if clusterColumns==1
    p=pdist(M.data',corr_type);
    ly=mylinkage(p,3);
    % [~,Yind,xx]=dendrogram(l,b);
    [~,Yind]=mydendrogram(ly,b);
    print -depsc dend.y.eps
else
    Yind=1:b;
end

M.data(nn)=NaN;
Mclust.data=M.data(Xind,Yind);

if diag==1
    [rowOrder, colOrder, newX] = diagonalize2d(Mclust.data);

    Xind=Xind(rowOrder);
    Yind=Yind(colOrder);

    Mclust.data=M.data(Xind,Yind);
    Mclust.rowlabels=M.rowlabels(Xind);
    Mclust.collabels=M.collabels(Yind);
else
    Mclust.data=M.data(Xind,Yind);
    Mclust.rowlabels=M.rowlabels(Xind);
    Mclust.collabels=M.collabels(Yind);
end

%%%
%%%
% min/max in colormap
%%%
%%%
% 
%  Mclust.data(Mclust.data>2)=2;
%  Mclust.data(Mclust.data<-2)=-2;


% don't plot if asked not to
if strcmp(varargin{end},'no')
    return
else
    heatmap2(Mclust.data,Mclust.collabels,Mclust.rowlabels,[],'colormap',...
    'hot',...
     'TickAngle',90,'ShowAllTicks', 1,'TickFontSize',10);
end

% figure out xind and yind
[~,Xind]=vlookup_list(Mclust.rowlabels,M.rowlabels,M.rowlabels);
[~,Yind]=vlookup_list(Mclust.collabels,M.collabels,M.collabels);
