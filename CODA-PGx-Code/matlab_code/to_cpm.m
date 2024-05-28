function cpm=to_cpm(x)

[ss,tt]=size(x);
cpm=nan(ss,tt);

for n=1:tt
    cpm(:,n)=x(:,n)./(sum(x(:,n)/1e6));
end
