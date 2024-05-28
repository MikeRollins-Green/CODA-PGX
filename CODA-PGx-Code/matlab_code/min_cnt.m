function pos=min_cnt(pos,min_cnt)

cnt=zeros(size(pos.collabels));
for nn=1:length(cnt)
cnt(nn)=abs(log2(median(pos.counts(pos.pool==1,nn))/median(pos.counts(pos.pool==2,nn))));
%cnt(nn)=log2(sum(pos.counts(pos.pool==1,nn))/sum(pos.counts(pos.pool==2,nn)));
end

pos=downSize(pos,cnt>=min_cnt);


