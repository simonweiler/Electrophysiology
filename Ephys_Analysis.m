%Analysis of slice electropysiology Margrie Lab started Spetember 2019
%% NOTES
%Data information is in excel sheet which can be called by Matlab 
%Nested functions
   % parseExperimentsXls_ephys
   % loadDataFile_wavesurfer
   % Pandora Toolbox 
   % spike event
   
%% Scripts starts here
%% Set directories and experimentator
experimentator= 'SW';
rdata_dir         = 'K:\SimonWeiler\RawData\Ephys_slice';%data directory of raw data;change accordingly
adata_dir         = 'K:\SimonWeiler\AnalyzedData\Ephys_slice';%data directory of extracted date;change accordingly
ExpXls            = 'C:\Users\Simon-localadmin\Documents\MargrieLab\Excel_exp_list\Experiment_list.xlsx';%directory where excel batch file is located;change accordingly
%% parse Experiments XLS database
batchopt          = parseExperimentsXls_ephys(ExpXls);%calls the nested function parseExperimentsXls_ephys and considers the user flag (1 or 0)
nummice           = length(batchopt.mouse);%length of experiments to be analyzed
 %% Go through each mouse and recordings iteratively
    
adder=1;%counting variable
 for i=1:nummice%for loop over experiments across days
     datapath=fullfile(rdata_dir, num2str(batchopt.mouse{i}), filesep);%directory and name of experiments (from excel sheet)
     cd(char(datapath));%go to directory
      for k=1:length(batchopt.exp_ids{i})%loop in bigger loop for each cell per experimental day
          if batchopt.exp_ids{i}(k)<10%for cells with id less then XX0010, e.g., XX0001-XX0009
                n_str = sprintf( '%04d', batchopt.exp_ids{i}(k));
          else
                n_str = sprintf( '%04d', batchopt.exp_ids{i}(k));%for cells with id mor then XX0010, e.g., XX0011-XX0099
          end
            fold_name=[experimentator n_str];%complete cell folder name such as SW0001 or MF0001
            exp_folder=fullfile(datapath,fold_name);%complete folder and directory
            cell_ab=[num2str(batchopt.mouse{i}) fold_name];
          % Check how many files are in the cell folder and cd to it
            list=dir([char(exp_folder) '\*.h5']);
            cd(exp_folder);
            
           %Check wether intrinsic properties or LED stimulation was
           %used 
            idx_led=find(contains({list(:).name},'LED')==1);
            intrlist={'IV','Passive','Rheobase','Passive','Sag'};
            idx_intr=find(contains({list(:).name},intrlist)==1);
            
            %Get LED stimuli of the cell
            if isempty(idx_led)==1
               disp('No LED evoked stimulus measured');
            else isempty(idx_led)==0
                for m=1:length(idx_led);  
                filename=[char(exp_folder) '\' list(idx_led(m)).name];
                %load file using loadDataFile_wavesurfer
                data = loadDataFile_wavesurfer(filename);
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
                %LED train 
                    if contains(stimuli_type,'train')==1 
                       delay=str2num(data.header.StimulusLibrary.Stimuli.element12.Delegate.Delay)*sr;
                       pulsedur=str2num(data.header.StimulusLibrary.Stimuli.element12.Delegate.PulseDuration)*sr;
                       endtime=data.header.StimulusLibrary.Stimuli.element12.Delegate.EndTime;
                       period=str2num(data.header.StimulusLibrary.Stimuli.element12.Delegate.Period);
                       duration=str2num(data.header.StimulusLibrary.Stimuli.element12.Delegate.Duration)*sr
                       
                       %stimvec=[firstamp:ampch:(pulsenr*ampch+(firstamp-ampch))];
                    else
                    end
                end
            end
            
            
            %Get the intrinsic properties of the cell
            if isempty(idx_intr)==1
                disp('No intrinsic properties measured');
            else isempty(idx_intr)==0
             for m=1:length(idx_intr);  
                filename=[char(exp_folder) '\' list(idx_intr(m)).name];
                %load file using loadDataFile_wavesurfer
                data = loadDataFile_wavesurfer(filename);
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
                 %IV display and extraction 
                    if contains(stimuli_type,'IV')==1 
                       delay=str2num(data.header.StimulusLibrary.Stimuli.element1.Delegate.Delay)*sr;
                       delaybp=str2num(data.header.StimulusLibrary.Stimuli.element1.Delegate.DelayBetweenPulses)*sr;
                       pulsenr=str2num(data.header.StimulusLibrary.Stimuli.element1.Delegate.PulseCount);
                       pulsedur=str2num(data.header.StimulusLibrary.Stimuli.element1.Delegate.PulseDuration)*sr;
                       firstamp=str2num(data.header.StimulusLibrary.Stimuli.element1.Delegate.FirstPulseAmplitude);
                       ampch=str2num(data.header.StimulusLibrary.Stimuli.element1.Delegate.AmplitudeChangePerPulse);
                       stimvec=[firstamp:ampch:(pulsenr*ampch+(firstamp-ampch))];
                       %Deconcatenate traces
                       startpoint=[delay/2:pulsedur+delaybp:tracel];
                       endpoint=[startpoint(2):pulsedur+delaybp:tracel];
                       for t=1:length(stimvec)
                       traces_deconc(:,t)=filt_traces(startpoint(t)+1:endpoint(t));
                       squarepulse(:,t)=[repmat(0,delay/2,1);repmat(stimvec(t),pulsedur,1);repmat(0,delay/2,1)];
                       end
                       %1 : Resting membrane potential;
                        for t=1:length(stimvec)
                        RM_tr(:,t)=traces_deconc(1:delay/2,t);
                        end
                        RMP=round(mean(RM_tr(:)),1);
                       %Plotting the IV and showing RMP
                       fig1=figure;set(fig1, 'Position', [200, 600, 200, 350]);set(gcf,'color','w');
                       subplot(2,1,1);
                       plot(traces_deconc,'Color','k','LineWidth',1);set(gca,'box','off');axis off;xlim([-260*srF length(traces_deconc)]);
                       hold on;text(-450*srF,RMP,[num2str(RMP),'mV'],'FontSize',9);
                       %Scale bar
                       scale_x= 200;
                       scale_y= 40;
                       ov_min=min(min(traces_deconc));
                       %scale barx
                       hold on;x1= (1250*srF)-(100*srF);x2=1250*srF;p1=plot([x1 x2],[ov_min-10 ov_min-10],'-','Color','k','LineWidth',1.5);
                       %scale bary
                       hold on;y2= (ov_min-10)+scale_y;y1=(ov_min-10);p2=plot([x2 x2],[y1 y2],'-','Color','k','LineWidth',1.5); 
                       title(cell_ab);
                       
                       subplot(2,1,2);
                       ov_min=min(min(squarepulse));
                       plot(squarepulse,'Color',[0.7 0.7 0.7],'LineWidth',1);set(gca,'box','off');axis off;xlim([-260*srF length(traces_deconc)]);
                       %Scale bar
                       scale_x= 200;
                       scale_y= 100;
                       %scale barx
                       hold on;x1= (1250*srF)-(100*srF);x2=1250*srF;p1=plot([x1 x2],[ov_min-10 ov_min-10],'-','Color','k','LineWidth',1.5);
                       %scale bary
                       hold on;y2= ov_min+scale_y;y1=ov_min;p1=plot([x2 x2],[y1 y2],'-','Color','k','LineWidth',1.5); 
                       %Spike event / frequency
                       for t=1:length(stimvec)
                       spike_event{:,t}=spike_times(traces_deconc(:,t),0.5);
                       spike_log(:,t)=~isempty(spike_event{:,t});
                       spikecount(:,t)=length(spike_event{:,t});
                       end
                       %Plotting spike frequency
                       fig2=figure;set(fig2, 'Position', [200, 800, 200, 200]);set(gcf,'color','w');
                       plot(stimvec,spikecount,'--bo');set(gca,'box','off');set(gcf,'color','w');xlabel('Current (pA)');ylabel('Spike Frequency (Hz)');
                       max_spikefr=max(spikecount);
                           
                 %Rheobase display and extraction     
                    elseif contains(stimuli_type,'Rheobase')==1;
                        delay=str2num(data.header.StimulusLibrary.Stimuli.element9.Delegate.Delay)*sr;
                        delaybp=str2num(data.header.StimulusLibrary.Stimuli.element9.Delegate.DelayBetweenPulses)*sr;
                        pulsenr=str2num(data.header.StimulusLibrary.Stimuli.element9.Delegate.PulseCount);
                        pulsedur=str2num(data.header.StimulusLibrary.Stimuli.element9.Delegate.PulseDuration)*sr;
                        firstamp=str2num(data.header.StimulusLibrary.Stimuli.element9.Delegate.FirstPulseAmplitude);
                        ampch=str2num(data.header.StimulusLibrary.Stimuli.element9.Delegate.AmplitudeChangePerPulse);
                        stimvec=[firstamp:ampch:(pulsenr*ampch+(firstamp-ampch))];
                        %Deconcatenate traces
                        startpoint=[delay/2:pulsedur+delaybp:tracel];
                        endpoint=[startpoint(2):pulsedur+delaybp:tracel];
                        for t=1:length(stimvec)
                        traces_deconc(:,t)=filt_traces(startpoint(t)+1:endpoint(t));
                        squarepulse(:,t)=[repmat(0,delay/2,1);repmat(stimvec(t),pulsedur,1);repmat(0,delay/2,1)];
                        end
                        for t=1:length(stimvec)
                        spike_event{:,t}=spike_times(traces_deconc(:,t),0.5);
                        spike_log(:,t)=~isempty(spike_event{:,t});
                        spikecount(:,t)=length(spike_event{:,t});
                        end
                        spike_idx=find(spike_log==1);
                        Rheobase=stimvec(spike_idx(1));
                        %Pandora parameters
                        trace_curr=traces_deconc(:,spike_idx(1));
                        dt = 1e-4;
                        dy = 1e-3;
                        props = struct('spike_finder', 2, 'threshold', 0);
                        traces_analysis=trace(trace_curr,dt, dy, 'Analysis', props);
                        alltrace_info = getProfileAllSpikes(traces_analysis);
                        parameters=alltrace_info.spikes_db.data; 
                        
                       %Plotting
                       fig4=figure;set(fig4, 'Position', [200, 600, 200, 300]);set(gcf,'color','w');
                       subplot(2,1,1);
                       plot(traces_deconc(:,spike_idx(1)),'Color','k','LineWidth',2);set(gca,'box','off');xlim([-260*srF length(traces_deconc)]);
                       hold on;text(-260*srF,round(parameters(3),1),[num2str(round(parameters(3),1)),'mV'],'FontSize',9);axis off;
                       %Scale bar
                       scale_x= 200;
                       scale_y= 40;
                       ov_min=min(min(traces_deconc(:,spike_idx(1))))-40;
                       %scale barx
                       hold on;x1= 1250*srF-(100*srF);x2=1250*srF;p1=plot([x1 x2],[ov_min-10 ov_min-10],'-','Color','k','LineWidth',1.5);
                       %scale bary
                       hold on;y2= (ov_min-10)+scale_y;y1=(ov_min-10);p2=plot([x2 x2],[y1 y2],'-','Color','k','LineWidth',1.5); 
                       title(cell_ab);
                       
                       subplot(2,1,2);
                       ov_min=min(min(squarepulse));
                       plot(squarepulse(:,spike_idx(1)),'Color',[0.7 0.7 0.7],'LineWidth',1);set(gca,'box','off');xlim([-260*srF length(traces_deconc)]);
                       axis off;
                       %Scale bar
                       scale_x= 200;
                       scale_y= 20;
                       %scale barx
                       hold on;x1= (1250*srF)-(100*srF);x2=1250*srF;p1=plot([x1 x2],[ov_min-30 ov_min-30],'-','Color','k','LineWidth',1.5);
                       %scale bary
                       hold on;y2= (ov_min-30)+scale_y;y1=(ov_min-30);p2=plot([x2 x2],[y1 y2],'-','Color','k','LineWidth',1.5); 
                       hold on;text(-260*srF,Rheobase,[num2str(round(Rheobase,1)),'pA'],'FontSize',9);axis off;
                    
                       
                       
                 %Passive display and extraction
                    elseif contains(stimuli_type,'Passive')==1; 
                           delay=str2num(data.header.StimulusLibrary.Stimuli.element6.Delegate.Delay)*sr;
                           delaybp=str2num(data.header.StimulusLibrary.Stimuli.element6.Delegate.DelayBetweenPulses)*sr;
                           pulsenr=str2num(data.header.StimulusLibrary.Stimuli.element6.Delegate.PulseCount);
                           pulsedur=str2num(data.header.StimulusLibrary.Stimuli.element6.Delegate.PulseDuration)*sr;
                           firstamp=str2num(data.header.StimulusLibrary.Stimuli.element6.Delegate.FirstPulseAmplitude);
                           ampch=str2num(data.header.StimulusLibrary.Stimuli.element6.Delegate.AmplitudeChangePerPulse);
                           stimvec=[firstamp:ampch:(pulsenr*ampch+(firstamp-ampch))];
                           %Deconcatenate traces
                           startpoint=[delay/2:pulsedur+delaybp:tracel];
                           endpoint=[startpoint(2):pulsedur+delaybp:tracel];
                           traces_deconc=[];
                           squarepulse=[];
                           for t=1:length(stimvec)
                           traces_deconc(:,t)=filt_traces(startpoint(t)+1:endpoint(t));
                           squarepulse(:,t)=[repmat(0,delay/2,1);repmat(stimvec(t),pulsedur,1);repmat(0,delay/2,1)];
                           end  
                           Vrest_IR=mean(traces_deconc(1:delay,:));
                           Vss_IR=mean(traces_deconc(2*delay:3*delay,:));
                           dV_IR=Vss_IR-Vrest_IR;
                           %Linear fit
                           P = polyfit(stimvec',dV_IR',1);
                           yfit = P(1)*stimvec+P(2);
                           %Plotting fit
                           fig3=figure;set(fig3, 'Position', [400, 800, 200, 200]);set(gcf,'color','w');
                           plot(stimvec,dV_IR,'ko');
                           hold on;ylabel('\DeltaVoltage (mV)');xlabel('Current (pA)');
                           plot(stimvec,yfit,'r-.');set(gca,'box','off');axis square;
                           %Input resistance
                           Rin=round(P(1)*1000,0);
                           hold on;
                           text(-30,max(dV_IR),[num2str(Rin),'M\Omega']);
                           %Sag display and extraction
                    elseif contains(stimuli_type,'Sag')==1;
                    %Need to implement then
                    else   
                        disp('Uncharted waters');
                    end
             end
            end
       end
  end