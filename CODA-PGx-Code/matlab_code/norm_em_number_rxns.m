function X=norm_em_number_rxns(pos2,min_bcg_counts,aa,bb)

% run through every drug sample, find proper bcg, and divide by bcg average
pos2=downSize(pos2,pos2.num_rxns>=aa&pos2.num_rxns<=bb);
oo=find(pos2.conc>0);

x=pos2;

for nn=1:length(oo)
    
    % find bcg
    bcg_off=strcmp(pos2.drug,'DMSO')&strcmp(pos2.cell_line,x.cell_line(oo(nn)))&x.pt_pool==pos2.pt_pool(oo(nn))&aa<=pos2.num_rxns(oo(nn))&bb>=pos2.num_rxns(oo(nn))&strcmp(x.plate,pos2.plate(oo(nn)));
    cc_bcg=pos2.data(:,bcg_off);
    bcg_mean=mean(cc_bcg,2);

    % replace outputs with nans where dmso is underpowered
    cc_bcg_counts=pos2.counts(:,bcg_off);
    bcg_mean_counts=mean(cc_bcg_counts,2);
    
    tt_dat=log2(pos2.data(:,oo(nn))./bcg_mean);
    tt_dat(bcg_mean_counts<min_bcg_counts)=nan;
    x.data(:,oo(nn))=tt_dat;
    
end

X=downSizeTo(x,oo,length(x.collabels));




