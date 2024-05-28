function [g,vaf]=sub2taste_genie_w_rna(ss)

% g=sub2gene(snv,in_domain,min_phast,min_cov)

% downsize to taste OR clinvar OR deogen>0.6 OR nonsense OR splice sites
% lin_off=ss.taste==1|ss.deogen>=0.6|ss.clinvar==1|ss.syn==3|ss.syn==4|ss.syn==5;

% lets ignore the frameshifts for now
%lin_off=ss.taste==1|ss.deogen>=0.6|ss.clinvar==1|ss.syn==3|ss.syn==5;
lin_off=ss.taste==1|ss.deogen>=0.6|ss.clinvar==1|ss.syn==3|ss.syn==4|ss.syn==5;
% lin_off=ss.taste>0.9|ss.deogen>=0.6|ss.clinvar==1|ss.syn==3|ss.syn==4|ss.syn==5;

% downsize on SNVs
ss=downSize(ss,lin_off);

% for purposes of ranking, stop (syn==3) is worse than indel (syn==4)
% since taking max below, convert 3 to 5 for this exercise only
% ss.syn(ss.syn==3)=5;

disp('done downsizing, now pulling genes');

% pull genes
[g.gm,bb,cc]=unique(ss.gm);
g.gene=ss.gene(bb);
g.collabels=ss.collabels;
g.data=sparse(length(g.gm),length(g.collabels));
g.altc=g.data;
g.refc=g.data;
g.rna_altc=g.data;
g.rna_refc=g.data;
g.aa=cell(length(g.gm),length(g.collabels));
vaf=g.data;
g.cov=g.data;

hh=1;
fprintf([num2str(max(cc)/10) ' dots to go' '\n']);
singles=0;
moret=0;

for ii=1:max(cc)
    
    gm_off=ii==cc;
    temp_dat=ss.data(gm_off,:)>0;
    temp_syn=ss.syn(gm_off);
    temp_aa=ss.aa(gm_off);
    temp_vaf=ss.data(gm_off,:);
    temp_cov=ss.cov(gm_off,:);
    temp_altc=ss.altc(gm_off,:);
    temp_refc=ss.refc(gm_off,:);
    temp_rna_altc=ss.rna_altc(gm_off,:);
    temp_rna_refc=ss.rna_refc(gm_off,:);
    
    for jj=1:length(ss.collabels)
        if sum(temp_dat(:,jj))>0
            [g.data(ii,jj),ioff]=max(temp_syn(temp_dat(:,jj)),[],1);
            xx=find(full(temp_dat(:,jj)));
            g.aa{ii,jj}=temp_aa{xx(ioff)};
            vaf(ii,jj)=temp_vaf(xx(ioff),jj);
            g.cov(ii,jj)=temp_cov(xx(ioff),jj);
            g.refc(ii,jj)=temp_refc(xx(ioff),jj);
            g.altc(ii,jj)=temp_altc(xx(ioff),jj);
            g.rna_refc(ii,jj)=temp_rna_refc(xx(ioff),jj);
            g.rna_altc(ii,jj)=temp_rna_altc(xx(ioff),jj);
                        
            if sum(temp_dat(:,jj))==1
                singles=singles+1;
            else
                moret=moret+1;
            end
            
            
        end
    end
    if hh==10
        fprintf('.');
        hh=1;
    end
    hh=hh+1;    
end

% back to conventional stop==3 notation
%g.data(g.data==5)=3;

g.vaf=vaf;

% disp(singles);
% disp(moret);




