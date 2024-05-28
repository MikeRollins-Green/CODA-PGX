function signi=calc_intra_sig2(posr,filt_expr,filt_guide_min)

% filt_guide_min = -1/0/1 (take all/ignore bad, keep no power/only keep
% validated sets)

% if dumping guides that don't cut
% 1= good guide, 0=no power, -1=bad guide (no cutting)
% of the ~180 genes, 25 become ineligible for readout b/c they have 3 or 4 bad guides

[guide,cutrate]=textread('C:\Users\Tomas\Dropbox\co\science\pipe\resources\guide_cutting_validation.txt','%s%n','delimiter','\t');
posr.guide_qc=vlookup(posr.rowlabels,guide,cutrate);
posr=downSize(posr,posr.guide_qc>=filt_guide_min);


% x dimensions
[~,tt]=size(posr.data);

% split into neg ctrls and drivers
ncr=downSize(posr,posr.pool==6);
posr=downSize(posr,posr.pool~=6);

% will need cc offsets later
[signi.gene,~,cc]=unique(posr.gene);

% initialize
signi.rowlabels=signi.gene;
signi.p=nan(max(cc),tt);
signi.fc=zeros(max(cc),tt);
signi.conc=posr.conc;
signi.cell_line=posr.cell_line;
signi.drug=posr.drug;
signi.collabels=posr.collabels;
signi.group=posr.group;
%signi.timepoint=posr.timepoint;

for jj=1:tt
    
    % for each gene compute a p-val
    for ii=1:max(cc)
        
        off=ii==cc;
        [signi.p(ii,jj),signi.fc(ii,jj)]=my_sig(posr.data(off,jj),ncr.data(:,jj));
    end    
end

% flip sign for GOF and plot
gof=textread('C:\Users\Tomas\Dropbox\co\science\queens_collab2\combined_data\gof_genes.txt','%s');
ii=ismember(signi.rowlabels,gof);
signi.fc(ii,:)=-signi.fc(ii,:);
signi.rowlabels=mark_gof(signi.rowlabels,ii);

if logical(filt_expr)
    % filter against mutations and expression
    load C:\Users\Tomas\Dropbox\co\science\pipe\lines_mut_expr\ll
    % expression level of "2" is a good threshold
    signi.filt=wes_rnaseq_filt(signi,ll,2);
end

signi.fc(signi.fc==-inf)=-3;
signi.fc(signi.fc==inf)=3;


%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%

function [p,fc]=my_sig(pools,ncr_pool)

pools=pools(isfinite(pools));

% init
p=nan;
fc=0;

if length(pools)>1
    
    % fc -> mean of drug/dmso b/w the two pools
    fc=mean(pools,1);
    
    % p - ranksum it    
    p=ranksum(pools,ncr_pool);
    
end



% old code for taking proportions to compute p-value
%
% for ii=1:ss
%
%     % fc -> mean of drug/dmso b/w the two pools
%     fc(ii)=mean([pool1(ii);pool2(ii)]);
%
%     % p_up
%     pt1=sum(ncr1>pool1(ii))/length(ncr1);
%     pt2=sum(ncr2>pool2(ii))/length(ncr2);
%     p_up=pt1*pt2;
%
%     % p_down
%     pt1=sum(ncr1<pool1(ii))/length(ncr1);
%     pt2=sum(ncr2<pool2(ii))/length(ncr2);
%     p_down=pt1*pt2;
%
%     if fc(ii)<0
%         p(ii)=p_down;
%     elseif fc(ii)>0
%         p(ii)=p_up;
%     end
%
% end




