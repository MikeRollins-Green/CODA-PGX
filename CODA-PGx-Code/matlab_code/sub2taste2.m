function [g,vaf]=sub2taste2(ss)

% g=sub2gene(snv,in_domain,min_phast,min_cov)

% downsize to taste OR syn==4 (deleterious indel/splice site)
lin_off=ss.taste==1|ss.syn==4;

% downsize on SNVs
ss=downSize(ss,lin_off);

% for purposes of ranking, stop (syn==3) is worse than indel (syn==4)
% since taking max below, convert 3 to 5 for this exercise only
ss.syn(ss.syn==3)=5;

disp('done downsizing, now pulling genes');

% pull genes
[g.gm,~,cc]=unique(ss.gm);
g.collabels=ss.collabels;
g.specimen=ss.specimen;
g.model=ss.model;
g.data=sparse(length(g.gm),length(g.collabels));
g.aa=cell(length(g.gm),length(g.collabels));
vaf=g.data;
g.cov=g.data;

hh=1;
fprintf([num2str(max(cc)/1000) ' dots to go' '\n']);
singles=0;
moret=0;

for ii=1:max(cc)
    
    gm_off=ii==cc;
    temp_dat=ss.data(gm_off,:)>0;
    temp_syn=ss.syn(gm_off);
    temp_aa=ss.aa(gm_off);
    temp_vaf=ss.data(gm_off,:);
    temp_cov=ss.cov(gm_off,:);
    
    for jj=1:length(ss.collabels)
        if sum(temp_dat(:,jj))>0
            [g.data(ii,jj),ioff]=max(temp_syn(temp_dat(:,jj)),[],1);
            xx=find(full(temp_dat(:,jj)));
            g.aa{ii,jj}=temp_aa{xx(ioff)};
            vaf(ii,jj)=temp_vaf(xx(ioff),jj);
            g.cov(ii,jj)=temp_cov(xx(ioff),jj);
                        
            if sum(temp_dat(:,jj))==1
                singles=singles+1;
            else
                moret=moret+1;
            end
            
            
        end
    end
    if hh==1000
        fprintf('.');
        hh=1;
    end
    hh=hh+1;    
end

% back to conventional stop==3 notation
g.data(g.data==5)=3;
g.vaf=vaf;

% disp(singles);
% disp(moret);




