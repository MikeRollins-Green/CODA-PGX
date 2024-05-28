function cc=string_em(pos)

cc={};

for n=1:length(pos.collabels)
    cc{n}=[pos.drug{n} '-' pos.line{n} '-' num2str(pos.ic(n)) '-' num2str(pos.day(n))];
end

cc=cc';
