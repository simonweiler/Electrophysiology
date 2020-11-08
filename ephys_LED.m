function [peak_n peak_p dpeak_n dpeak_p sub_traces] = ephys_LED(filename, data, clamp,amplifier)

fc=7;
base_start=1;
resp_win=0.05%50 ms in s

sr=data.header.StimulationSampleRate;
if str2num(filename(end-4))==0
temp=str2num(filename(end-3));
if amplifier==1
traces=data.(['sweep_000',num2str(temp)]).analogScans(:,1);
else amplifier==2
traces=data.(['sweep_000',num2str(temp)]).analogScans(:,2);
end
else
temp=str2num([filename(end-4) filename(end-3)]);
if amplifier==1
traces=data.(['sweep_00',num2str(temp)]).analogScans(:,1);
else amplifier==2
traces=data.(['sweep_00',num2str(temp)]).analogScans(:,2);
end
end

srF=sr/1000;               
       %Filter traces
        cutoff      = 1000;     % Hz (use 500 Hz for mini event / amplitude detection and 1000Hz for max currents. Chen & Regehr 2000)
        order       = 4;        % filter order ('pole'). (use 4 pole for minis and max current. Chen & Regehr 2000)
        type        = 'Butter';
        filt_traces = lowpassfilt(traces, order, cutoff, sr, type);
        tracel=length(traces);
        stimuli_type=data.header.DataFileBaseName;

  %LED TRAIN/SQUARE
 
    if contains(stimuli_type,'train')==1 || contains(stimuli_type,'square')==1
         bs_std=[];   
      bs_traces=[];  
      neg_peak=[];
      neg_fail=[];
      pos_peak=[];
      pos_fail=[];
       delay=str2num(data.header.StimulusLibrary.Stimuli.element12.Delegate.Delay)*sr;
       pulsedur=str2num(data.header.StimulusLibrary.Stimuli.element12.Delegate.PulseDuration)*sr;
       endtime=data.header.StimulusLibrary.Stimuli.element12.Delegate.EndTime*sr;
       period=str2num(data.header.StimulusLibrary.Stimuli.element12.Delegate.Period)*sr;
       duration=str2num(data.header.StimulusLibrary.Stimuli.element12.Delegate.Duration)*sr
       freq=duration/period;                                               
       %deconcatenate traces
       ste=[1:period:endtime];
       for k=1:freq
      traces_clip=  filt_traces(ste(k):ste(k)+0.5*sr,:); 
      bs=traces_clip(delay-pulsedur:delay,:);%first 100 ms baseline trace
      bs_std(k)=std(bs);%std of baseline trace
      bs_traces(:,k)=traces_clip-mean(traces_clip(base_start:delay,:));%subtract baseline
      
       %Negative peak
          neg_peak(k)=min(bs_traces(delay:delay+0.05*sr,k));
          neg_fail(k)=neg_peak(k)<fc*bs_std(k)*(-1);
          tr_d=bs_traces(delay:delay+0.05*sr,k);
          delay_peak_n(k)=find(tr_d==neg_peak(k));
          %Positive peak
          pos_peak(k)=max(bs_traces(delay:delay+0.05*sr,k));
          pos_fail(k)=pos_peak(k)>fc*bs_std(k); 
          delay_peak_p(k)=find(tr_d==pos_peak(k));
          %Replace values with 0 where 3std is not applicable
           if neg_fail(k)==1;
           peak_train_n(k)= neg_peak(k);
           dpeak_n(k)=delay_peak_n(k);
           else neg_fail(k)==0;
           peak_train_n(k)= 0;
           dpeak_n(k)=0;
           end
           if pos_fail(k)==1;
           peak_train_p(k)= pos_peak(k);
           dpeak_p(k)=delay_peak_p(k);
           else pos_fail(k)==0;
           peak_train_p(k)= 0;
           dpeak_p(k)=0;
           end 
       end  
       
             peak_n=peak_train_n;
             peak_p=peak_train_p;
             sub_traces=filt_traces-mean(filt_traces(base_start:delay,:));
             sub_traces=sub_traces(:,:);
   %LADDER  
     elseif contains(stimuli_type,'ladder')==1;
      bs_std=[];   
      bs_traces=[];  
      neg_peak=[];
      neg_fail=[];
      pos_peak=[];
      pos_fail=[];
      dpeak_n=[];
      dpeak_p=[];
      sub_traces=[];
      delay=str2num(data.header.StimulusLibrary.Stimuli.element13.Delegate.Delay)*sr;
      pulsedur=str2num(data.header.StimulusLibrary.Stimuli.element13.Delegate.PulseDuration)*sr;
      endtime=data.header.StimulusLibrary.Stimuli.element13.Delegate.EndTime*sr;
      delay_b_p=str2num(data.header.StimulusLibrary.Stimuli.element13.Delegate.DelayBetweenPulses)*sr;
      pulse_count=str2num(data.header.StimulusLibrary.Stimuli.element13.Delegate.PulseCount);
      %deconcatenate traces
      ind_traces=reshape(filt_traces(1:endtime,:),[endtime/pulse_count pulse_count]);
      %run through each trace
      for k=1:pulse_count
          traces_clip=ind_traces(:,k);
          bs=traces_clip(delay-10*pulsedur:delay,:);%first 100 ms baseline trace
          bs_std(k)=std(bs);%std of baseline trace
          bs_traces(:,k)=traces_clip-mean(traces_clip(base_start:delay,:));%subtract baseline
          
          %Negative peak
          neg_peak(k)=min(bs_traces(delay:end,k));
          neg_fail(k)=neg_peak(k)<fc*bs_std(k)*(-1); 
          dpeak_n(k)=NaN;
          %Positive peak
          pos_peak(k)=max(bs_traces(delay:end,k));
          pos_fail(k)=pos_peak(k)>fc*bs_std(k); 
          dpeak_p(k)=NaN;
          %Replace values with 0 where 3std is not applicable
           if neg_fail(k)==1;
           peak_ladder_n(k)= neg_peak(k);
           else neg_fail(k)==0;
           peak_ladder_n(k)= 0;
           end
           if pos_fail(k)==1;
           peak_ladder_p(k)= pos_peak(k);
           else pos_fail(k)==0;
           peak_ladder_p(k)= 0;
           end       
         
                        
      end
           
             peak_n=peak_ladder_n;
             peak_p=peak_ladder_p;
             sub_traces=filt_traces-mean(filt_traces(base_start:delay,:));
             sub_traces=sub_traces(:,:);
            
     
    elseif contains(stimuli_type,'high')==1 
         bs_std=[];   
      bs_traces=[];  
      neg_peak=[];
      neg_fail=[];
      pos_peak=[];
      pos_fail=[];
       delay=str2num(data.header.StimulusLibrary.Stimuli.element18.Delegate.Delay)*sr;
       pulsedur=str2num(data.header.StimulusLibrary.Stimuli.element18.Delegate.PulseDuration)*sr;
       endtime=data.header.StimulusLibrary.Stimuli.element18.Delegate.EndTime*sr;
       if contains(stimuli_type,'high2')==1
       period=str2num(data.header.StimulusLibrary.Stimuli.element19.Delegate.Period)*sr;
       else
       period=str2num(data.header.StimulusLibrary.Stimuli.element18.Delegate.Period)*sr;
       end
       duration=str2num(data.header.StimulusLibrary.Stimuli.element18.Delegate.Duration)*sr;
       freq=duration/period;                                               
       %deconcatenate traces
       ste=[1:period:endtime];
       for k=1:freq
      traces_clip=  filt_traces(ste(k):ste(k)+0.5*sr,:); 
      bs=traces_clip(delay-pulsedur:delay,:);%first 100 ms baseline trace
      bs_std(k)=std(bs);%std of baseline trace
      bs_traces(:,k)=traces_clip-mean(traces_clip(base_start:delay,:));%subtract baseline
      
       %Negative peak
          neg_peak(k)=min(bs_traces(delay:delay+0.05*sr,k));
          neg_fail(k)=neg_peak(k)<fc*bs_std(k)*(-1);
          tr_d=bs_traces(delay:delay+0.05*sr,k);
          delay_peak_n(k)=find(tr_d==neg_peak(k));
          %Positive peak
          pos_peak(k)=max(bs_traces(delay:delay+0.05*sr,k));
          pos_fail(k)=pos_peak(k)>fc*bs_std(k); 
          delay_peak_p(k)=find(tr_d==pos_peak(k));
          %Replace values with 0 where 3std is not applicable
           if neg_fail(k)==1;
           peak_train_n(k)= neg_peak(k);
           dpeak_n(k)=delay_peak_n(k);
           else neg_fail(k)==0;
           peak_train_n(k)= 0;
           dpeak_n(k)=0;
           end
           if pos_fail(k)==1;
           peak_train_p(k)= pos_peak(k);
           dpeak_p(k)=delay_peak_p(k);
           else pos_fail(k)==0;
           peak_train_p(k)= 0;
           dpeak_p(k)=0;
           end 
       end  
       
             peak_n=peak_train_n;
             peak_p=peak_train_p;
             sub_traces=filt_traces-mean(filt_traces(base_start:delay,:));
             sub_traces=sub_traces(:,:);
           
    else
             peak_n=[];
             peak_p=[];
             sub_traces=[];
     end                 
 %CC or VC  

 
end