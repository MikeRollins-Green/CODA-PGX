function pos=mk_pos(my_dir)

drop_samples_lt_mean_reads_guide=-1; % these get removed later in the dmso min cov math

[pos.original_sort,pos.gwz_id,pos.collabels,pos.operator,pos.cell_line,pos.pt_pool,pos.drug,pos.conc,pos.group]=textread([my_dir '\collabels.3.txt'],'%n%s%s%s%s%n%s%n%n','delimiter','\t');
[pos.rowlabels,pos.gene,pos.pool]=textread([my_dir '\rowlabels.txt'],'%s%s%n','delimiter','\t');
pos.counts=load([my_dir '\matrix.tab']);
pos.counts(isnan(pos.counts))=0;
pos.cpm=to_cpm(pos.counts);
pos.pool(isnan(pos.pool))=6;
pos.data=quantilenorm(pos.cpm);

% drop the crap pools
incl_these=sum(pos.counts)./length(pos.rowlabels)>=drop_samples_lt_mean_reads_guide;
fprintf(['dropping ' num2str(sum(~incl_these)) ' samples since dont have at least ' num2str(drop_samples_lt_mean_reads_guide) ' reads per guide' '\n']);
pos=downSize(pos,incl_these);
