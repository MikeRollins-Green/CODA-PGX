function X=norm_em_poc(nc,min_bcg_counts)

% run through every drug sample, find proper bcg, and divide by bcg average
x=nc;
oo=find(nc.conc>0);

for nn=1:length(oo)
    
    % find bcg
    bcg_off=strcmp(nc.drug,'DMSO')&strcmp(nc.cell_line,x.cell_line(oo(nn)))&x.pt_pool==nc.pt_pool(oo(nn))&x.timepoint==nc.timepoint(oo(nn))&strcmp(x.plate,nc.plate(oo(nn)));
    cc_bcg=nc.data(:,bcg_off);
    bcg_mean=mean(cc_bcg,2);

    % replace outputs with nans where dmso is underpowered
    cc_bcg_counts=nc.counts(:,bcg_off);
    bcg_mean_counts=mean(cc_bcg_counts,2);
    
    tt_dat=log2(nc.data(:,oo(nn))./bcg_mean);

    % keep stag2 and brca1
    ook=find(strcmp('BRCA1',nc.gene)|strcmp('STAG2',nc.gene));
    k_dat=tt_dat(ook,:);
    
    tt_dat(bcg_mean_counts<min_bcg_counts)=nan;
    tt_dat(ook,:)=k_dat;
    tt_dat(bcg_mean_counts<50)=nan;
    
    
    x.data(:,oo(nn))=tt_dat;
    
end

X=downSizeTo(x,oo,length(x.collabels));




