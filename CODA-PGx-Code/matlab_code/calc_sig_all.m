function signi=calc_sig_all(posr,ncr)

% compute a fold-change and significance matrix for each drug x conc x cell
% line, compare all concentrations for a drug (could have more pool-values going into
% ranksum test)

% initialize
signi.gene=posr.gene;
signi.rowlabels=signi.gene;

signi.p=[];
signi.fc=[];

signi.drug=[];
signi.collabels=[];
signi.cell_line=[];

% for dumping the redundant half of samples
uniq_tag=[];

% for each sample
[~,tt]=size(posr.data);

for nn=1:tt
    
    % oo=find(posr.conc==posr.conc(nn)&strcmp(posr.drug,posr.drug(nn))&strcmp(posr.cell_line,posr.cell_line(nn)));
    oo=find(strcmp(posr.drug,posr.drug(nn)));
    temp_pools=posr.pt_pool(oo);
    
    % only use samples where you have both pools read out
    if sort(unique(temp_pools))==[1;2]
        
        % the data
        [temp_p,temp_fc]=my_sig(posr.data(:,oo),ncr.data(:,oo));
        signi.p=[signi.p,temp_p];
        signi.fc=[signi.fc,temp_fc];
        
        % the labels
        signi.cell_line=[signi.cell_line;posr.cell_line(nn)];
        signi.drug=[signi.drug;posr.drug(nn)];
        signi.collabels=[signi.collabels;posr.drug(nn)];
        
        uniq_tag=[uniq_tag; {[char(posr.drug(nn))]}];
        
    end
end

[~,oo]=unique(uniq_tag);
signi=downSizeTo(signi,oo,length(signi.collabels));

% flip sign for GOF and plot
gof=textread('C:\Users\Tomas\Dropbox\co\science\queens_collab2\combined_data\gof_genes.txt','%s');
ii=ismember(signi.rowlabels,gof);
signi.fc(ii,:)=-signi.fc(ii,:);
signi.rowlabels=mark_gof(signi.rowlabels,ii);
% filter against mutations and expression
load C:\Users\Tomas\Dropbox\co\science\pipe\lines_mut_expr\ll
signi.filt=wes_rnaseq_filt(signi,ll,2);





function [p,fc]=my_sig(pos_data,ncr_data)

% merge nc pools and do a single ranksum test
% probably a major violation of stats laws... check with
% statistician on this

ncr_data=zscore(ncr_data);
%pos_data=zscore(pos_data)

ncr_pool=ncr_data(:);

% init
[ss,~]=size(pos_data);

p=ones(ss,1);
fc=zeros(ss,1);

for ii=1:ss
    
    % fc -> mean of drug/dmso b/w the two pools
    fc(ii)=nanmean(pos_data(ii,2));
    
    % p - ranksum it
    
    if sum(isnan(pos_data(ii,:)),2)>0
        p(ii)=nan;
    else
        p(ii)=ranksum(pos_data(ii,:)',ncr_pool);
    end
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




