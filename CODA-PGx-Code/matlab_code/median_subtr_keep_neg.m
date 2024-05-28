function d=median_subtr_keep_neg(d)

d(isnan(d))=0;

% median center
m=median(d,2);
for n=1:length(m)
    d(n,:)=d(n,:)-m(n);
%     x=d(n,:)<0;
%     d(n,x)=0;    
end
