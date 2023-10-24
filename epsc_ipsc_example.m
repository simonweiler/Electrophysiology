function epsc_ipsc_example(ex_trace,in_trace,sample_rate,show_time,start_time,time_stop,yborders,titlab)
fig4=figure;set(fig4, 'Position', [200, 800, 300, 300]);set(gcf,'color','w');
p1=plot([1:1:size(show_time,2)]*(1/sample_rate),ex_trace(show_time,1),'Color','r','linestyle', '-','LineWidth',1.2);hold on; ylim(yborders);p1.Color(4)=1;
p1=plot([1:1:size(show_time,2)]*(1/sample_rate),in_trace(show_time,1)+50,'Color','b','linestyle', '-','LineWidth',1.2);hold on; ylim(yborders);p1.Color(4)=1;

hold on;line([(start_time-show_time(1))*(1/sample_rate), ((start_time-show_time(1))+time_stop)*(1/sample_rate)], [yborders(2)-100 yborders(2)-100], 'color', 'c', 'linestyle', '-','LineWidth',0.5);hold on;
hold on;box off;axis off;hold on;title(titlab);set(gca,'FontSize',12);
hold on;line([9000*(1/sample_rate), 9500*(1/sample_rate)], [yborders(1)+50 yborders(1)+50], 'color', 'k', 'linestyle', '-','LineWidth',1);hold on;
 hold on;line([9000*(1/sample_rate), 9000*(1/sample_rate)], [yborders(1)+50 yborders(1)+150], 'color', 'k', 'linestyle', '-','LineWidth',1);hold on;
end