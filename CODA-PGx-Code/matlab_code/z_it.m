function x=z_it(in,dims)

x=zeros(length(in),dims);

for nn=1:dims
    x(:,nn)=in;
end
