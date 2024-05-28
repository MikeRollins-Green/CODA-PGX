function x=downSizeSpec(x,oo)

x_backup=x;
x=downSize(x,oo);
x.rowlabels=x_backup.rowlabels;
x.gene=x_backup.gene;
x.pool=x_backup.pool;

x.data=x_backup.data(:,oo);
x.cpm=x_backup.cpm(:,oo);
x.counts=x_backup.counts(:,oo);


