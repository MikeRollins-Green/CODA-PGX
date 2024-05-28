function X=norm_em_dmsoXdmso(pos2,min_bcg_counts)

% run through every drug sample, find proper bcg, and divide by bcg average
x=downSize(pos2,strcmp(pos2.drug,'DMSO'));
oo=length(x.collabels);

for nn=1:oo
    
    % find bcg
    bcg_off=logical(ones(oo,1));
    bcg_off(nn)=0;
    
    cc_bcg=pos2.data(:,bcg_off);
    bcg_mean=mean(cc_bcg,2);

    % replace outputs with nans where dmso is underpowered
    cc_bcg_counts=pos2.counts(:,bcg_off);
    bcg_mean_counts=mean(cc_bcg_counts,2);
    
    tt_dat=log2(x.data(:,nn)./bcg_mean);
    tt_dat(bcg_mean_counts<min_bcg_counts)=nan;
    x.data(:,nn)=tt_dat;
    
end

X=x;




