function view_raw44(my_gene,cnum,posr,pos,min_guide_val)

% nn= column offset in posr
% min guide val -1 = all, 0=skip -1 which are bad guides, 1 = only take
% guides validated

nn=find(strcmp(cnum,posr.collabels));

if isempty(nn)
    disp('your drug is not in this dataset');
    return
end


[guide,cutrate]=textread('C:\Users\TomasBabak\Dropbox (personal)\co\science\pipe\resources\guide_cutting_validation.txt','%s%n','delimiter','\t');
guide_qc=vlookup(posr.rowlabels,guide,cutrate);

posr=downSize(posr,guide_qc>=min_guide_val);
pos=downSize(pos,guide_qc>=min_guide_val);

abs_edg=0:1:5000;

ii=strcmp(my_gene,posr.gene);

if sum(ii)==0
    disp('cant find your gene -> make sure it exists in posr.gene');
else    

    % the data - ratio neg ctrls, ratio guide of BRCA1
    ncrat=posr.data(strcmp(posr.gene,'NEGCTRL'),nn);
    sgrat=posr.data(ii,nn);
    
    aa=max(abs([ncrat;sgrat]));
%    ratio_edg=-10:0.01:10;
    ratio_edg=-aa*1.5:aa/500:aa*1.5;

    % raw dat
    dmso_off=strcmp(pos.drug,'DMSO')&strcmp(pos.cell_line,posr.cell_line(nn));
    sg_off=strcmp(posr.collabels(nn),pos.collabels);    
    
    % neg controls - my_sample, DMSO
    ncsg=pos.counts(strcmp(posr.gene,'NEGCTRL'),sg_off);
    ncdmso=pos.counts(strcmp(posr.gene,'NEGCTRL'),dmso_off);
    
    % BRCA1 - my_sample, DMSO
    sgsg=pos.counts(ii,sg_off);
    sgdmso=pos.counts(ii,dmso_off);
    
    % absolute cpm pool 2, dmso
    subplot(3,1,3)
    hold off
    aa=histc(mean(ncdmso,2),abs_edg);
    bar(abs_edg,aa,1);
    ylim([0 1]);
    hold on
    aa=histc(mean(sgdmso,2),abs_edg);
    bar(abs_edg,aa,4,'r');
    ylim([0 1]);
    title(['DMSO neg ctrl counts (blue) ||  ' my_gene ' counts in (red)']);
    xlabel('counts')
    set(gca,'YTickLabel',[]);
     
    % absolute cpm pool 2, dmso
    subplot(3,1,2)
    hold off
    aa=histc(ncsg,abs_edg);
    bar(abs_edg,aa,1);
    ylim([0 1]);
    hold on
    aa=histc(sgsg,abs_edg);
    bar(abs_edg,aa,4,'r');
    ylim([0 1]);
    title([cnum ' treated neg ctrl counts (blue) || ' my_gene ' counts in (red)']);
    xlabel('counts')
    set(gca,'YTickLabel',[]);
    
    % absolute cpm pool 2, dmso
    subplot(3,1,1)
    hold off
    aa=histc(ncrat,ratio_edg);
    bar(ratio_edg,aa,1,'b');
    ylim([0 1]);
    hold on
    aa=histc(sgrat,ratio_edg);
    bar(ratio_edg,aa,4,'r');
    ylim([0 1]);
    title([cnum ' treated neg ctrl drug/dmso (blue) || ' my_gene ' ratios (red)']);
    xlabel('log2(ratio)')
    set(gca,'YTickLabel',[]);
  
       
end



