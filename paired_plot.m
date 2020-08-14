function paired_plot(data,test_u)
%plot indivudal date points with lines connecting them and plot median+ do
%test
fig3= figure;set(fig3, 'Name', 'Paired comp');set(fig3, 'Position', [200, 300, 200, 200]);set(gcf,'color','w');
hold on;
for i=1:length(data)
     p1=plot([1,2],[data(:,1),data(:,2)],'color',[0.5 0.5 0.5]);    
end
hold on;pS=plotSpread([data(:,1),data(:,2)],'categoryIdx',[ones(1,length(data(:,1)))' ones(1,length(data(:,2)))'*2],...
    'categoryMarkers',{'o','o'},'categoryColors',{'k','r'});hold on;
hold on;plot([1,2],[nanmedian(data(:,1)),nanmedian(data(:,2))],'k','LineWidth',3);
box off;set(gca,'FontSize',10);
if test_u==0
[p1]=signrank(data(:,1) ,data(:,2));title([' p=' num2str(p1) ', n=' num2str(length(data))],'FontWeight','Normal');
else
 [p1]=ranksum(data(:,1) ,data(:,2));title([' p=' num2str(p1) ', n=' num2str(length(data))],'FontWeight','Normal');
end
end