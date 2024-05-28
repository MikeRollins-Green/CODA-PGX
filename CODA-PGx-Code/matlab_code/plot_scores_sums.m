function plot_scores_sums(pos)

pp=downSize(pos,strcmp(pos.plate,'60')&strcmp(pos.drug,'DMSO'))
pp.collabels=pp.num_rxns
subplot(2,2,1)
mc=corr_all_no_plot(pp)
plot(mc.collabels,sum(mc.data),'o')
xlabel('number of reactions pooled')
ylabel('sum of correlations with all reactions')
grid on
title('DMSO 60mm')

pp=downSize(pos,strcmp(pos.plate,'60')&strcmp(pos.drug,'Olaparib'))
pp.collabels=pp.num_rxns
subplot(2,2,2)
mc=corr_all_no_plot(pp)
plot(mc.collabels,sum(mc.data),'o')
xlabel('number of reactions pooled')
ylabel('sum of correlations with all reactions')
grid on
title('Olap 60mm')

pp=downSize(pos,strcmp(pos.plate,'T75')&strcmp(pos.drug,'DMSO'))
pp.collabels=pp.num_rxns
subplot(2,2,3)
mc=corr_all_no_plot(pp)
plot(mc.collabels,sum(mc.data),'o')
xlabel('number of reactions pooled')
ylabel('sum of correlations with all reactions')
grid on
title('DMSO T75')

pp=downSize(pos,strcmp(pos.plate,'T75')&strcmp(pos.drug,'Olaparib'))
pp.collabels=pp.num_rxns
subplot(2,2,4)
mc=corr_all_no_plot(pp)
plot(mc.collabels,sum(mc.data),'o')
xlabel('number of reactions pooled')
ylabel('sum of correlations with all reactions')
grid on
title('Olap T75')