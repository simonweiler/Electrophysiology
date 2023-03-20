%readout spikecounts/frequency

function [maxcount_spk] = spike_counts(Ephys,cell_log)
temp=[];
temp=find(cell_log==1);
for i=1:length(temp)
    if  isempty(Ephys(temp(i)).IV)==0
    maxcount_spk(i)=max(Ephys(temp(i)).IV.spikecount);
    else
    maxcount_spk(i)=NaN;
    end
end

end 