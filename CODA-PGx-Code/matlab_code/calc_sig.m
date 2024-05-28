function signi=calc_sig(posr,ncr)

% compute a fold-change and significance matrix for each drug x conc x cell
% line

% initialize
signi.gene=posr.gene;
signi.rowlabels=signi.gene;

signi.p=[];
signi.fc=[];

signi.conc=[];
signi.cell_line=[];
signi.drug=[];
signi.collabels=[];

% for dumping the redundant samples
uniq_tag=[];

% for each sample
[~,tt]=size(posr.data);

for nn=1:tt
    
    oo=find(posr.conc==posr.conc(nn)&strcmp(posr.drug,posr.drug(nn))&strcmp(posr.cell_line,posr.cell_line(nn)));
    temp_pools=posr.pt_pool(oo);
    
    % only use samples where you have at least two pools reporting
    if length(unique(temp_pools))>1
        
        % the data
        [temp_p,temp_fc]=my_sig(posr.data(:,oo(1)),posr.data(:,oo(2)),ncr.data(:,oo(1)),ncr.data(:,oo(2)));
        signi.p=[signi.p,temp_p];
        signi.fc=[signi.fc,temp_fc];
        
        % the labels
        signi.conc=[signi.conc;posr.conc(nn)];
        signi.drug=[signi.drug;posr.drug(nn)];
        signi.cell_line=[signi.cell_line;posr.cell_line(nn)];
        signi.collabels=[signi.collabels;posr.collabels(nn)];
        
        uniq_tag=[uniq_tag; {[char(num2str(posr.conc(nn))) char(posr.drug(nn)) char(posr.cell_line(nn))]}];
        
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

signi.fc(signi.fc==-inf)=-3;
signi.fc(signi.fc==inf)=3;





function [p,fc]=my_sig(pool1,pool2,ncr1,ncr2)

% merge nc pools and do a single ranksum test
% probably a major violation of stats laws... check with
% statistician on this
ncr_pool=[ncr1;ncr2];

% init
p=ones(size(pool1));
fc=zeros(size(pool1));

ss=length(p);

for ii=1:ss
    
    % fc -> mean of drug/dmso b/w the two pools
    fc(ii)=mean([pool1(ii);pool2(ii)]);
    
    % p - ranksum it
    
    if isnan(pool1(ii))|isnan(pool2(ii))
        p(ii)=nan;
    else
        p(ii)=ranksum([pool1(ii);pool2(ii)],ncr_pool);
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




