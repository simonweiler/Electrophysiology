function plot_spike_freq(data,idx)
fig5= figure;set(fig5, 'Name', 'Plot3');set(fig5, 'Position', [200, 100, 400, 400]);set(gcf,'color','w');
for i=1:length(find(idx==1));
    temp=find(idx==1);
   if isempty(data(temp(i)).IV)==0
    plot(data(temp(i)).IV.stimvec,data(temp(i)).IV.spikecount,'--ok'); 
    
    hold on;
    else  
      plot(NaN);
    end
end

ylabel('Spike frequency (Hz');xlabel('Current (pA)');box off
end
