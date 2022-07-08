function [bs_traces pulse_specs] = ephys_spikeligth(filename, data)
base_start=1;
resp_win=0.05;%50 ms in s

sr=data.header.StimulationSampleRate;
if str2num(filename(end-4))==0
temp=str2num(filename(end-3));
traces=data.(['sweep_000',num2str(temp)]).analogScans(:,1);
else
temp=str2num([filename(end-4) filename(end-3)]);
traces=data.(['sweep_00',num2str(temp)]).analogScans(:,1);

end

srF=sr/1000;               
       %Filter traces
        cutoff      = 1000;     % Hz (use 500 Hz for mini event / amplitude detection and 1000Hz for max currents. Chen & Regehr 2000)
        order       = 4;        % filter order ('pole'). (use 4 pole for minis and max current. Chen & Regehr 2000)
        type        = 'Butter';
        filt_traces = lowpassfilt(traces, order, cutoff, sr, type);
        tracel=length(traces);
        stimuli_type=data.header.DataFileBaseName;

        delay=str2num(data.header.StimulusLibrary.Stimuli.element21.Delegate.Delay)*sr;
       pulsedur=str2num(data.header.StimulusLibrary.Stimuli.element21.Delegate.PulseDuration)*sr;
       endtime=data.header.StimulusLibrary.Stimuli.element21.Delegate.EndTime*sr;
       period=str2num(data.header.StimulusLibrary.Stimuli.element21.Delegate.Period)*sr;
       duration=str2num(data.header.StimulusLibrary.Stimuli.element21.Delegate.Duration)*sr;
       amplit=str2num(data.header.StimulusLibrary.Stimuli.element21.Delegate.Amplitude);
       freq=duration/period; 

if  ismember(stimuli_type,"LED_before_rheo_1")
       pdurLED=str2num(data.header.StimulusLibrary.Stimuli.element25.Delegate.PulseDuration)*sr;
       ampLED=str2num(data.header.StimulusLibrary.Stimuli.element25.Delegate.Amplitude);
        delayLED=str2num(data.header.StimulusLibrary.Stimuli.element25.Delegate.Delay)*sr;
        type_stim=1;
elseif ismember(stimuli_type,"LED_before_rheo_2")
       pdurLED=str2num(data.header.StimulusLibrary.Stimuli.element26.Delegate.PulseDuration)*sr;
       ampLED=str2num(data.header.StimulusLibrary.Stimuli.element26.Delegate.Amplitude);
       delayLED=str2num(data.header.StimulusLibrary.Stimuli.element26.Delegate.Delay)*sr;
       type_stim=2;
elseif ismember(stimuli_type,"LED_during_rheo_1")
       pdurLED=str2num(data.header.StimulusLibrary.Stimuli.element23.Delegate.PulseDuration)*sr;
       ampLED=str2num(data.header.StimulusLibrary.Stimuli.element23.Delegate.Amplitude);
       delayLED=str2num(data.header.StimulusLibrary.Stimuli.element23.Delegate.Delay)*sr;
       type_stim=3;
else ismember(stimuli_type,"LED_during_rheo_2")
    pdurLED=str2num(data.header.StimulusLibrary.Stimuli.element24.Delegate.PulseDuration)*sr;
       ampLED=str2num(data.header.StimulusLibrary.Stimuli.element24.Delegate.Amplitude);
       delayLED=str2num(data.header.StimulusLibrary.Stimuli.element24.Delegate.Delay)*sr;
       type_stim=4;
end

bs_traces=[];traces_clip=[];bs=[];bs_std=[];
       ste=[1:period:endtime];
       for k=1:freq
      traces_clip=  filt_traces(ste(k):ste(k)+1.5*sr,:); 
       bs=traces_clip(1:delay-0.15*sr,:);%first 100 ms baseline trace
       bs_std(k)=std(bs);%std of baseline trace
      bs_traces(:,k)=traces_clip-mean(traces_clip(1:delay-0.15*sr,:));%subtract baseline
       end
       
pulse_specs.stimtype=type_stim;
pulse_specs.amp_current=amplit;
pulse_specs.durLED=pdurLED;
pulse_specs.delayLED=delayLED;


end
