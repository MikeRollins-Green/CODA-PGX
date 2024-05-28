function [out,keep,my_mem]=merge_concentrations(pos)

out=pos;

oo=find(pos.conc>0);
my_mem=zeros(length(pos.collabels));

for nn=1:length(oo)

    % update counts
    oo2=find(strcmp(pos.drug,pos.drug(oo(nn)))&pos.pt_pool==pos.pt_pool(oo(nn))&strcmp(pos.cell_line,pos.cell_line(oo(nn))));

    %%
    for jj=1:length(oo2)
         out.counts(:,oo2(jj))=sum(pos.counts(:,oo2),2);
    end
    
    % ultimately just want one entry per drug
    my_mem(oo(nn),oo2)=1;
    
end

% keep
[ss,~]=size(my_mem);
keep=[];

for ii=1:ss
    
    toff=find(my_mem(ii,:)>0);

    if isempty(toff)
        keep=[keep; ii];
    else
        keep=[keep; min(toff)];
    end
end
keep=unique(keep);

out=downSizeTo(out,keep,length(out.collabels));
out.conc(out.conc>0)=10;
    
    
    
    
    
    



