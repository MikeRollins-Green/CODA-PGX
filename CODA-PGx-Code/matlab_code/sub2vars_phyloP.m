function g=sub2vars_phyloP(ss,my_min_syn,my_max_gnom,my_min_revel,my_min_vest,my_max_fathmm,my_min_mpc,my_min_taste,my_min_clinvar,my_min_deogen,my_min_clinpred,my_min_in_domain,my_min_phast,my_min_phyloP)

% g=sub2gene(snv,in_domain,min_phast,min_cov)

% downsize to our params, OR syn==4 (deleterious indel) OR syn==5 (change in
% canonical splice site)
lin_off=ss.gnom<=my_max_gnom&...
    ss.revel>=my_min_revel&...
    ss.vest>=my_min_vest&...
    ss.fathmm<=my_max_fathmm&...
    ss.mpc>=my_min_mpc&...
    ss.taste>=my_min_taste&...
    ss.phyloP>=my_min_phyloP&...
    ismember(ss.clinvar,my_min_clinvar)&...
    ss.deogen>=my_min_deogen&...
    ss.clinpred>=my_min_clinpred&...
    ismember(ss.syn,my_min_syn)&...
    ss.phast>=my_min_phast&...
    ismember(ss.in_domain,my_min_in_domain);

% include deletions
% lin_off=lin_off|ss.syn==3|ss.syn==4|ss.syn==5;

% downsize on SNVs
ss=downSize(ss,lin_off);

% for purposes of ranking, stop (syn==3) is worse than indel (syn==4)
% since taking max below, convert 3 to 6 for this exercise only
ss.syn(ss.syn==3)=6;

disp('done downsizing, now pulling genes');

% pull genes
[g.gm,~,cc]=unique(ss.gm);
g.collabels=ss.collabels;
g.tissue=ss.tissue;
g.id=ss.id;
g.data=sparse(length(g.gm),length(ss.id));

hh=1;
fprintf([num2str(max(cc)/1000) ' dots to go' '\n']);

for ii=1:max(cc)
    
    gm_off=ii==cc;
    temp_dat=ss.data(gm_off,:)>0;
    temp_syn=ss.syn(gm_off);
    
    for jj=1:length(ss.id)
        if sum(temp_dat(:,jj))>0
            %disp(jj);
            g.data(ii,jj)=max(temp_syn(temp_dat(:,jj)),[],1);
        end
    end
    if hh==1000
        fprintf('.');
        hh=1;
    end
    hh=hh+1;    
end

% back to conventional stop==3 notation
g.data(g.data==6)=3;





