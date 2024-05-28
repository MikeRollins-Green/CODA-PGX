function plot_volcanoes_q_pdf(signi,output_dir)

signi.p=signi.q;

for n=1:length(signi.collabels)
    
    expr=signi.filt(:,n)>0;
    gof=signi.filt(:,n)==2;
    
    hold off
    plot(signi.fc(expr,n),-log10(signi.p(expr,n)),'o','markersize',10,'linewidth',1.5,'color',[18/255 179/255 101/255])
    grid on
    hold on
    plot(signi.fc(gof,n),-log10(signi.p(gof,n)),'*','markersize',10,'linewidth',2,'color',[14/255 66/255 60/255])
    xmin_ss=max(abs(signi.fc(expr,n)))*1.2;
    %     if xmin_ss>3
    %         xmin_ss=3;
    %     end
    text_off=(signi.p(:,n)<0.5&abs(signi.fc(:,n))>0.01&expr)|gof;
    text(signi.fc(text_off,n),-log10(signi.p(text_off,n)),signi.rowlabels(text_off),'fontsize',6,'color',[14/255 66/255 60/255])
    title(signi.collabels(n),'fontsize',10);
    xlabel('log2(Fold-Change)','fontsize',10)
    ylabel('-log10(FDR)','fontsize',10)
    
    if isnan(xmin_ss)
        continue
    else
        xlim([-xmin_ss xmin_ss])
    end
    grid on;
    set(gca,'fontsize',10)
    %print('-djpeg',[output_dir '/' signi.collabels{n} '.jpg'],'-r120')
    print('-dpdf',[output_dir '/' signi.collabels{n} '.pdf'],'-r120')
    fprintf([num2str(sum(signi.fc(n,:)==-inf)) '\n']);
end
