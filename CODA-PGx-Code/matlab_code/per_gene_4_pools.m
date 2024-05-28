function pos3=per_gene_4_pools(pos2)

off1=find(pos2.pool==1&strcmp(pos2.gene,'NEGCTRL')~=1);
off2=find(pos2.pool==2&strcmp(pos2.gene,'NEGCTRL')~=1);
off3=find(pos2.pool==3&strcmp(pos2.gene,'NEGCTRL')~=1);
off4=find(pos2.pool==4&strcmp(pos2.gene,'NEGCTRL')~=1);

[pos3.rowlabels,aa,~]=intersect(pos2.rowlabels(off1),pos2.rowlabels(off2));
pos3.gene=pos2.gene(off1(aa));

% seed the data matrices - 4x the size of input
pos3.data=pos2.data(off1(aa),:);
pos3.cpm=pos2.cpm(off1(aa),:);
pos3.counts=pos2.counts(off1(aa),:);

[~,~,bb]=intersect(pos2.gene(off1),pos2.gene(off2));
pos3.data(:,pos2.pt_pool==2)=pos2.data(off2(bb),pos2.pt_pool==2);
pos3.cpm(:,pos2.pt_pool==2)=pos2.cpm(off2(bb),pos2.pt_pool==2);
pos3.counts(:,pos2.pt_pool==2)=pos2.counts(off2(bb),pos2.pt_pool==2);

pos3.collabels=pos2.collabels;
pos3.operator=pos2.operator;
pos3.cell_line=pos2.cell_line;
pos3.pt_pool=pos2.pt_pool;
pos3.drug=pos2.drug;
pos3.conc=pos2.conc;
pos3.timepoint=pos2.timepoint;
pos3.plate=pos2.plate;


