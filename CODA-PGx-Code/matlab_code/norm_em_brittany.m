function X=norm_em_brittany(nc,min_bcg_counts)

% run through every drug sample, find proper bcg, and divide by bcg average
x=nc;
oo=find(nc.conc>0);

for nn=1:length(oo)
    
    % find bcg
    bcg_off=strcmp(nc.drug,'DMSO')&strcmp(nc.cell_line,x.cell_line(oo(nn)))&x.pt_pool==nc.pt_pool(oo(nn));
    cc_bcg=nc.data(:,bcg_off);
    bcg_mean=mean(cc_bcg,2);
    
    fid=fopen(['data_out/' x.cell_line{oo(nn)} '-' num2str(nc.pt_pool(oo(nn))) '.tab'],'w');
    
    for iii=1:length(bcg_mean)
        fprintf(fid,[num2str(bcg_mean(iii)) '\n']);
    end
    fclose(fid);
    
    
    % replace outputs with nans where dmso is underpowered
    cc_bcg_counts=nc.counts(:,bcg_off);
    bcg_mean_counts=mean(cc_bcg_counts,2);
    
    tt_dat=log2(nc.data(:,oo(nn))./bcg_mean);
    tt_dat(bcg_mean_counts<min_bcg_counts)=nan;
    x.data(:,oo(nn))=tt_dat;
    
end

X=downSizeTo(x,oo,length(x.collabels));




