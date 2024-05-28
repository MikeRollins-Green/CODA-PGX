function [out1,out2]=prune_posr(posr,ncr)

% make sure both pools exist for each cell line, conc, drug combo


[ss,tt]=size(posr.data);
keep=logical(zeros(tt,1));

for nn=1:tt
    
    oo=posr.conc==posr.conc(nn)&strcmp(posr.drug,posr.drug(nn))&strcmp(posr.cell_line,posr.cell_line(nn));
    temp_pools=posr.pt_pool(oo);
    if sort(temp_pools)==[1;2]
        keep(nn)=1;
    end
    
end

out1=downSize(posr,keep);
out2=downSize(ncr,keep);


