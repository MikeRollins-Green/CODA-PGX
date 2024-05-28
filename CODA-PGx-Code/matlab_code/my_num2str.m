function yy=my_num2str(xx)

yy={};
for nn=1:length(xx)
    yy{nn}=num2str(xx(nn));
end
