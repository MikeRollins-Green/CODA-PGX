function cycle_pools(pos)

[~,tt]=size(pos.data);
% pos.cpm(pos.cpm==0)=1;

for n=1:tt
    
    %
    my_boxplot(pos.cpm(pos.pool==0,n),pos.cpm(pos.pool==1,n),pos.cpm(pos.pool==2,n),pos.cpm(pos.pool==3,n),pos.cpm(pos.pool==4,n),pos.cpm(pos.pool==5,n),pos.cpm(pos.pool==6,n));
    
    xticks([1 2 3 4 5 6])
    xticklabels({'pool1','pool2','pool3','pool4','yulei_pool','pool1/2_neg_ctrl'})
    title(pos.collabels{n});
    set(gca,'yscale','log')
    ylabel('CPM')
    grid on
    ylim([1 1e4])
    waitforbuttonpress;
    
end
