function my_breakdown(pos)

pos.data(isnan(pos.data))=0;
[~,xx,~,~,exp]=pca(pos.data');
figure(1)
bar(exp)
xlabel('Principal Components')
ylabel('Variance Explained')
title('Principal Component Analysis')
% assocation with variables

figure(2)

plot_it(xx,pos.pt_pool,'Pool',1)
plot_it(xx,pos.operator,'Operator',3)
plot_it(xx,pos.conc,'Concentation',5)
plot_it(xx,pos.drug,'Drug',7)
plot_it(xx,pos.cell_line,'Cell Line',9)
plot_it(xx,pos.pcr_cycles,'PCR Cycles',11)


function plot_it(xx,my_var,plot_title,nn)

[~,tt]=size(xx);
rr=zeros(tt,1);
pp=ones(tt,1);

subplot(6,2,nn)
% operator
[~,~,cc]=unique(my_var);
for n=1:tt
[rr(n),pp(n)]=corr(xx(:,n),cc);
end
bar(rr)
ylabel({'Var of PC explained by ';plot_title})
xlabel('Principal Component');
subplot(6,2,nn+1)
bar(-log10(pp))
ylabel('-log10(p) plot on left')
xlabel('Principal Component');

