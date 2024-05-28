function yy=my_num2cell(xx)

yy={};

for nn=1:length(xx)
    yy{nn}={num2str(xx(nn))};
end


