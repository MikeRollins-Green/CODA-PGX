function [J,p]=analyze_em(drugs,clin,x,my_genes,my_drug,my_cancer,max_trt_line)

%my_cancer=strcmp(drugs.cancer_type,'BRCA')|strcmp(drugs.cancer_type,'OV');

% pull patients with drug, cancer type, and max treatment line
off=strcmp(drugs.drug,my_drug)&strcmp(drugs.cancer_type,my_cancer)&drugs.treatment_order<=max_trt_line;
%off=strcmp(drugs.drug,my_drug)&drugs.treatment_order<=max_trt_line&my_cancer;
pts=drugs.sample(off);

% subset pts into mut/non-mut
[~,aa]=intersect(x.sample,pts);
x=downSizeTo(x,aa,length(x.sample));
[~,aa]=intersect(x.gene,my_genes);
x=downSizeTo(x,aa,length(x.gene));

J=nan;
p=nan;
muts_off=sum(x.data,1)>0;
if sum(muts_off)>0
      
    % km
    [J,p]=km_plot_two_samples_emerald(clin,x.sample(muts_off),x.sample(~muts_off));
    legend('muts','censored','wt');
    
end
