function paired_subplot(par1,par2,test_u)

str={'V_{min}(mV)','V_{peak}(mV)','V_{init}(mV)','V_{thresh}(mV)', 'Vslope_{max} (\DeltamV/\Deltams)','V_{half} (mV)','Spike_{amplitude} (mV)',...
    'AHP_{max}(mV)', 'Spike init (ms)','Spike_{rise} (ms)', 'Spike_{fall} (ms)','Spike_{base width} (ms)','Spike_{half width} (ms)'};
xlab=({'CPN','NtsR1'});

fig3= figure;set(fig3, 'Name', 'Paired comp');set(fig3, 'Position', [200, 300, 800, 800]);set(gcf,'color','w');
for i=1:size(par1,1)
data=[];
data=[par1(i,:)' par2(i,:)'];
subplot(4,4,i);hold on;
for k=1:length(data)
    p1=[];
    p1=plot([1,2],[data(:,1),data(:,2)],'color',[0.5 0.5 0.5]);    
    
end

hold on;pS=plotSpread([data(:,1),data(:,2)],'categoryIdx',[ones(1,length(data(:,1)))' ones(1,length(data(:,2)))'*2],...
    'categoryMarkers',{'o','o'},'categoryColors',{'k','r'});hold on;
hold on;plot([1,2],[nanmedian(data(:,1)),nanmedian(data(:,2))],'k','LineWidth',3);
box off;set(gca,'FontSize',10);
hold on;ylabel(str{i});xticklabels(xlab)
if test_u==0
[p1]=signrank(data(:,1) ,data(:,2));title([' p=' num2str(p1) ', n=' num2str(length(data))],'FontWeight','Normal');
else
 [p1]=ranksum(data(:,1) ,data(:,2));title([' p=' num2str(p1) ', n=' num2str(length(data))],'FontWeight','Normal');
end
end
end
