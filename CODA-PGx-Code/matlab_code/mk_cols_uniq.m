function xx=mk_cols_uniq(xx)

for n=1:length(xx)
    xx{n}=[xx{n} '-' num2str(n)];
end
