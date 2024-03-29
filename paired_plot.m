function stats=paired_plot(data,test_u,cl)
%plot indivudal date points with lines connecting them and plot median+ do
%test
fig3= figure;set(fig3, 'Name', 'Paired comp');set(fig3, 'Position', [200, 300, 200, 200]);set(gcf,'color','w');
hold on;
for i=1:length(data)
     pl=plot([1,2],[data(:,1),data(:,2)],'color',[0.8 0.8 0.8]);    
end

hold on;pS=plotSpread([data(:,1),data(:,2)],'categoryIdx',[ones(1,length(data(:,1)))' ones(1,length(data(:,2)))'*2],...
    'categoryMarkers',{'o','o'},'categoryColors',cl);hold on;
%hold on;plot([1,2],[nanmedian(data(:,1)),nanmedian(data(:,2))],'k','LineWidth',3);
hold on;plot([1,2],[nanmean(data(:,1)),nanmean(data(:,2))],'k','LineWidth',3);
box off;set(gca,'FontSize',10);
if test_u==0
[p1]=signrank(data(:,1) ,data(:,2));p1=round(p1,3);
title([' p=' num2str(p1) ', n=' num2str(length(data))],'FontWeight','Normal');
stats=p1;
elseif test_u==1
   [p1]=signrank(data(:,1) ,data(:,2),'tail','left');p1=round(p1,3);
title([' p=' num2str(p1) ', n=' num2str(length(data))],'FontWeight','Normal');
stats=p1;
elseif test_u==2
       [p1]=signrank(data(:,1) ,data(:,2),'tail','right');p1=round(p1,3);
title([' p=' num2str(p1) ', n=' num2str(length(data))],'FontWeight','Normal');
stats=p1;
else
 [p1]=ranksum(data(:,1) ,data(:,2));p1=round(p1,3);title([' p=' num2str(p1) ', n=' num2str(length(data))],'FontWeight','Normal');
 stats=p1;
end
end