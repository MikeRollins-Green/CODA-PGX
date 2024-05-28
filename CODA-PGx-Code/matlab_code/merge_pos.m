function pos=merge_pos(pos1,pos2)

% check to make sure collabels are unique
if ~isempty(intersect(pos1.collabels,pos2.collabels))
    disp('your collabels are not unique');
    return;
end


pos.collabels=[pos1.collabels;pos2.collabels];
pos.gwz=[pos1.gwz;pos2.gwz];
pos.operator=[pos1.operator;pos2.operator];
pos.cell_line=[pos1.cell_line;pos2.cell_line];
pos.pt_pool=[pos1.pt_pool;pos2.pt_pool];
pos.drug=[pos1.drug;pos2.drug];
pos.conc=[pos1.conc;pos2.conc];
pos.pcr_cycles=[pos1.pcr_cycles;pos2.pcr_cycles];

pos.rowlabels=unique([pos1.rowlabels;pos2.rowlabels]);

[~,aa1,bb1]=intersect(pos.rowlabels,pos1.rowlabels);
[~,aa2,bb2]=intersect(pos.rowlabels,pos2.rowlabels);

pos.gene=pos.rowlabels;

pos.gene(aa1)=pos1.gene(bb1);
pos.gene(aa2)=pos2.gene(bb2);

pos.pool=zeros(size(pos.rowlabels));

pos.pool(aa1)=pos1.pool(bb1);
pos.pool(aa2)=pos2.pool(bb2);

[~,tt1]=size(pos1.counts);

pos.counts=zeros(length(pos.pool),length(pos.collabels));
pos.counts(aa1,1:tt1)=pos1.counts(bb1,:);
pos.counts(aa2,tt1+1:end)=pos2.counts(bb2,:);
