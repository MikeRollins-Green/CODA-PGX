function signi2survival_inverse_fc(signi,clin,drugs,gof,lof,gg,min_fc,max_p_signi,minNpatients,max_treat_line,max_p_km,look,min_abs_J)

% signi2survival(signi,slin,drugs,gof,lof,min_fc,max_p_signi,minNpatients,max_treat_line,max_p_km)
% min_fc (e.g. 2 = less than -2 in log2 space)
% max_p_signi (e.g. 0.1 -> will only input genes which have dropout of p<0.1)
% minNpatients (min number of patients treated with drug and recorded in TCGA)
% max_treat_line (e.g. 1 = only first line therapy considered)
% max_p_km (only print significant hits to this p-value)

[~,tt]=size(signi.p);

% for each drug, check and print all the combos

for nn=1:tt
    
    my_drug=signi.drug(nn);
    
    % skip if drug isnt in tcga
    if isempty(intersect(my_drug,look.mydrug))
        %fprintf([my_drug{1} ' is not in TCGA' '\n']);
        continue
    end
    
    tcga_drug=vlookup_list(my_drug,look.mydrug,look.tcga_drug);
    
    % what cancer types have the patients (check up to 3)
    [tc,nc,tc2,nc2,tc3,nc3]=print_drug_usage_no_print(tcga_drug,drugs);
    
    % genes
    my_genes=signi.gene(signi.fc(:,nn)>min_fc&signi.p(:,nn)<=max_p_signi);
    
    if nc>=minNpatients
        run_em(drugs,clin,lof,my_genes,tcga_drug,tc,max_treat_line,max_p_km,'LOF',min_abs_J);
        run_em(drugs,clin,gof,my_genes,tcga_drug,tc,max_treat_line,max_p_km,'GOF',min_abs_J);
        run_em(drugs,clin,gg,my_genes,tcga_drug,tc,max_treat_line,max_p_km,'MUT_ONLY',min_abs_J);
    end
    
    if nc2>=minNpatients
        
        tc=tc2;
        run_em(drugs,clin,lof,my_genes,tcga_drug,tc,max_treat_line,max_p_km,'LOF',min_abs_J);
        run_em(drugs,clin,gof,my_genes,tcga_drug,tc,max_treat_line,max_p_km,'GOF',min_abs_J);
        run_em(drugs,clin,gg,my_genes,tcga_drug,tc,max_treat_line,max_p_km,'MUT_ONLY',min_abs_J);
    end
    
    if nc3>=minNpatients
        
        tc=tc3;
        run_em(drugs,clin,lof,my_genes,tcga_drug,tc,max_treat_line,max_p_km,'LOF',min_abs_J);
        run_em(drugs,clin,gof,my_genes,tcga_drug,tc,max_treat_line,max_p_km,'GOF',min_abs_J);
        run_em(drugs,clin,gg,my_genes,tcga_drug,tc,max_treat_line,max_p_km,'MUT_ONLY',min_abs_J);
    end
    
    
    
end


function run_em(drugs,clin,lof,my_genes,my_drug,tc,max_treat_line,max_p_km,my_ptype,min_abs_J)

% all genes
[J,p]=analyze_em_no_km(drugs,clin,lof,my_genes,my_drug,tc,max_treat_line);
if p<=max_p_km && abs(J)>min_abs_J
    for nn=1:length(my_genes)-1
        fprintf([my_genes{nn} ',']);
    end
    fprintf([my_genes{end} '\t' char(my_ptype) '\t' char(my_drug) '\t' char(tc) '\t' num2str(max_treat_line) '\t' num2str(p) '\t' num2str(J) '\n']);
end

% the other genes individually
for nn=1:length(my_genes)-1
    
    [J,p]=analyze_em_no_km(drugs,clin,lof,my_genes{nn},my_drug,tc,max_treat_line);
    if p<=max_p_km && abs(J)>min_abs_J
        fprintf([my_genes{nn} '\t' char(my_ptype) '\t' char(my_drug) '\t' char(tc) '\t' num2str(max_treat_line) '\t' num2str(p) '\t' num2str(J) '\n']);
    end
    
end





