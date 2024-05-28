function x=map_labels(n1,n2,n3,n4)

x={};

for nn=1:length(n1)
    x{nn}=[num2str(n4(nn)) '-' num2str(n1(nn)) '-' n2{nn} '-' n3{nn}];
end

