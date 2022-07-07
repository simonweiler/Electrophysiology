function ephys_spikeligth(filename, data)
sr=data.header.StimulationSampleRate;
temp=str2num(filename(end-3));
traces=data.(['sweep_000',num2str(temp)]).analogScans(:,1);
srF=sr/1000;               
       %Filter traces
        cutoff      = 1000;     % Hz (use 500 Hz for mini event / amplitude detection and 1000Hz for max currents. Chen & Regehr 2000)
        order       = 4;        % filter order ('pole'). (use 4 pole for minis and max current. Chen & Regehr 2000)
        type        = 'Butter';
        filt_traces = lowpassfilt(traces, order, cutoff, sr, type);
        tracel=length(traces);
        stimuli_type=data.header.DataFileBaseName;

        delay=str2num(data.header.StimulusLibrary.Stimuli.element12.Delegate.Delay)*sr;
       pulsedur=str2num(data.header.StimulusLibrary.Stimuli.element12.Delegate.PulseDuration)*sr;
       endtime=data.header.StimulusLibrary.Stimuli.element12.Delegate.EndTime*sr;
       period=str2num(data.header.StimulusLibrary.Stimuli.element12.Delegate.Period)*sr;
       duration=str2num(data.header.StimulusLibrary.Stimuli.element12.Delegate.Duration)*sr
       freq=duration/period;   

end
