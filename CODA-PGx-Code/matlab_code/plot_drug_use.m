function plot_drug_use(my_drug,drugs)

off=strcmp(drugs.drug,my_drug);
ct=unique(drugs.cancer_type);

ss=zeros(length(ct),1);
for n=1:length(ss)
    ss(n)=sum(strcmp(drugs.cancer_type(off),ct(n)));
end

[~,ii]=sort(ss,'descend');
bar(ss(ii));
set(gca,'xtick',[1:length(ct)]);
set(gca,'xticklabel',ct(ii),'xticklabelrotation',90)


