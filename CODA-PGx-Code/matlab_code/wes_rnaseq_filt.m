function x=wes_rnaseq_filt(signi,ll,min_expr)

% for each cell line x gene, mark x as: 0=not expressed, 1=all good (wt
% expressed), 2=gof mutated line

[ss,tt]=size(signi.p);

x=ones(ss,tt);
% labels
[~,aa,bb]=intersect(signi.gene,ll.gene);

for jj=1:tt
    
    % find line
    oo=find(strcmp(signi.cell_line(jj),ll.line));
    
    % not expr
    not_expr=ll.expr_data(:,oo)<min_expr;
    
    % gof mut
    gof_mut=ll.mut(:,oo)==1&ll.gof;
    
    my_dat=ones(length(ll.gene),1);
    my_dat(gof_mut)=2;
    my_dat(not_expr)=0;    
    
    x(aa,jj)=my_dat(bb);
end

    
