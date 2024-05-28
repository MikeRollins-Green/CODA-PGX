function pos2=duplicate_pool3(pos2)

oo=find(pos2.pt_pool==3);

if ~isempty(oo)
    pos2.collabels=[pos2.collabels;pos2.collabels(oo)];
    pos2.operator=[pos2.operator;pos2.operator(oo)];
    pos2.cell_line=[pos2.cell_line;pos2.cell_line(oo)];
    pos2.pt_pool(oo)=1;
    pos2.pt_pool=[pos2.pt_pool;ones(length(oo),1).*2];
    pos2.drug=[pos2.drug;pos2.drug(oo)];
    pos2.conc=[pos2.conc;pos2.conc(oo)];
    pos2.data=[pos2.data,pos2.data(:,oo)];
    pos2.counts=[pos2.counts,pos2.counts(:,oo)];
    pos2.cpm=[pos2.cpm,pos2.cpm(:,oo)];
end

    
    