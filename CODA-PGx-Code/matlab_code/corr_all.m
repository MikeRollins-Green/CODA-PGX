function mc=corr_all(pos)

%pos.data(pos.data<-5)=-5;
%pos.data(pos.data>5)=5;
%pos.data(isnan(pos.data))=0;

x.data=corr(pos.data);
%x.data=corr(quantilenorm(pos.cpm));
x.collabels=pos.collabels;
x.rowlabels=pos.collabels;

close all; mc=myCluster(x,-1,0,1,1,'correlation');
colormap hot;
axis square
colorbar;

