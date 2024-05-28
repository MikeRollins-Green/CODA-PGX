function X=norm_em3(pos2,min_bcg_counts)

% run through every drug sample, find proper bcg, and divide by bcg average
x=pos2;
oo=find(pos2.conc>0);

for nn=1:length(oo)
    
    % find bcg
    %bcg_off=strcmp(pos2.drug,'DMSO')&strcmp(pos2.cell_line,x.cell_line(oo(nn)))&x.pt_pool==pos2.pt_pool(oo(nn))&x.timepoint==pos2.timepoint(oo(nn))&strcmp(x.plate,pos2.plate(oo(nn)));
    bcg_off=x.group==pos2.group(oo(nn))&strcmp(pos2.drug,'DMSO')&strcmp(pos2.cell_line,x.cell_line(oo(nn)))&x.pt_pool==pos2.pt_pool(oo(nn));
    cc_bcg=pos2.data(:,bcg_off);
    bcg_mean=mean(cc_bcg,2);

    % replace outputs with nans where dmso is underpowered
    cc_bcg_counts=pos2.counts(:,bcg_off);
    bcg_mean_counts=mean(cc_bcg_counts,2);
    
    tt_dat=log2(pos2.data(:,oo(nn))./bcg_mean);
    tt_dat(bcg_mean_counts<min_bcg_counts)=nan;
    x.data(:,oo(nn))=tt_dat;
    
end
X=x;
%X=downSizeTo(x,oo,length(x.collabels));




