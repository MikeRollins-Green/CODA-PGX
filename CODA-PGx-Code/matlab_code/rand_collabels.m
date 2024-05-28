function xx=rand_collabels(in)

[~,ii]=sort(rand(length(in),1));

oo=1:length(in);
oo=oo(ii);

xx={};

for n=1:length(oo)
    xx{n}=['q' num2str(oo(n))];
end

        
    