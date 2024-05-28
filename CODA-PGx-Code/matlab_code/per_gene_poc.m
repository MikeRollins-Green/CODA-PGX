function pos3=per_gene_poc(pos2)

off1=find(pos2.pool==1&strcmp(pos2.gene,'NEGCTRL')~=1);
off2=find(pos2.pool==2&strcmp(pos2.gene,'NEGCTRL')~=1);
off1=find(pos2.pool==1&strcmp(pos2.gene,'NEGCTRL')~=1);
off2=find(pos2.pool==2&strcmp(pos2.gene,'NEGCTRL')~=1);


[pos3.rowlabels,aa,bb]=intersect(pos2.gene(off1),pos2.gene(off2));
pos3.gene=pos2.gene(off1(aa));
% copy from pos2
pos3.data=pos2.data(off1(aa),:);
pos3.cpm=pos2.cpm(off1(aa),:);
pos3.counts=pos2.counts(off1(aa),:);
%% update the segments where actually want pool 2
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


