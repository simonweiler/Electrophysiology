function epsc_ipsc_example(ex_trace,in_trace,sample_rate,show_time,start_time,time_stop,yborders,titlab)

%plots example EPSC and IPSC for a given cell

fig4=figure;set(fig4, 'Position', [200, 800, 300, 300]);set(gcf,'color','w');
%plot traces
p1=plot([1:1:size(show_time,2)]*(1/sample_rate),ex_trace(show_time,1),'Color','r','linestyle', '-','LineWidth',1.8);hold on; ylim(yborders);p1.Color(4)=1;
p1=plot([1:1:size(show_time,2)]*(1/sample_rate),in_trace(show_time,1)+50,'Color','b','linestyle', '-','LineWidth',1.8);hold on; ylim(yborders);p1.Color(4)=1;

hold on;line([(start_time-show_time(1))*(1/sample_rate), ((start_time-show_time(1))+time_stop)*(1/sample_rate)], [yborders(2)-100 yborders(2)-100], 'color', 'c', 'linestyle', '-','LineWidth',1);hold on;
hold on;box off;axis off;hold on;title(titlab,'FontWeight','normal');set(gca,'FontSize',12);
%scale bar
hold on;line([9000*(1/sample_rate), 9500*(1/sample_rate)], [yborders(1)+50 yborders(1)+50], 'color', 'k', 'linestyle', '-','LineWidth',1);hold on;
 hold on;line([9000*(1/sample_rate), 9000*(1/sample_rate)], [yborders(1)+50 yborders(1)+150], 'color', 'k', 'linestyle', '-','LineWidth',1);hold on;

 hold on;text(0-1600*(1/sample_rate),max(ex_trace(show_time,1)-10),'-70 mV','FontSize',9,'Color','r');
 hold on;text(0-1250*(1/sample_rate),max(ex_trace(show_time,1)+60),'0 mV','FontSize',9,'Color','b');

 hold on;text(5000*(1/sample_rate),min(ex_trace(show_time,1)-50),'EPSC','FontSize',9,'Color','r');
 hold on;text(5000*(1/sample_rate),max(in_trace(show_time,1)+50),'IPSC','FontSize',9,'Color','b');
end