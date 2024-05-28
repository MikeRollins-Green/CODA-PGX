function xx=mark_gof(in,ii)

% put + after GOF, - after LOF

xx=in;

iii=find(ii);

for n=1:length(iii)
    xx{iii(n)}=[in{iii(n)} '+'];
end

ii=~ii;
iii=find(ii);

for n=1:length(iii)
    xx{iii(n)}=[in{iii(n)} '-'];
end
