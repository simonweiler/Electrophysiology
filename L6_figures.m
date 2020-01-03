
%L6 plotting scripts 
%% Load L6 structure 
str_L6    = 'K:\SimonWeiler\AnalyzedData\Ephys_slice';
folder_list = uipickfiles('FilterSpec',str_L6);
load(char(folder_list));
srF=20;
%% Intrinsic properties
%% Read out overall max and min for figure 
step=2;
for i=1:size(Ephys,2)
    if isempty(Ephys(i).IV)==0
     min_tr(:,i)=min(min(Ephys(i).IV.traces(:,1:step:end))); 
     max_tr(:,i)=max(max(Ephys(i).IV.traces(:,1:step:end))); 
    else
     min_tr(:,i)=NaN;  
    end
end
ov_min=nanmin(min_tr);
ov_max=nanmax(max_tr);
%% Plotting the IV and showing RMP
fig1=figure;set(fig1, 'Position', [200, -200, 1200, 1200]);set(gcf,'color','w');
for i=1:size(Ephys,2)
hold on;subplot(7,6,i);
if isempty(Ephys(i).IV)==0
%Plot only every second trace atm for display
plot(Ephys(i).IV.traces(:,1:2:end),'Color','k','LineWidth',1);set(gca,'box','off');axis off;
%xlim([-260*srF length(traces_deconc)]);
hold on;ylim([ov_min-10 ov_max]);
hold on;text(-550*srF,Ephys(i).IV.RMP,[num2str(Ephys(i).IV.RMP),'mV'],'FontSize',9);
 %Scale bar
 scale_x= 200;
 scale_y= 40;
 %scale barx
 hold on;x1= (1250*srF)-(100*srF);x2=1250*srF;p1=plot([x1 x2],[ov_min-10 ov_min-10],'-','Color','k','LineWidth',1.5);
 %scale bary
 hold on;y2= (ov_min-10)+scale_y;y1=(ov_min-10);p2=plot([x2 x2],[y1 y2],'-','Color','k','LineWidth',1.5); 
title([num2str(Ephys(i).patching_date) Ephys(i).celID]);
else
plot(1,1);set(gca,'box','off');axis off;
end
end
%% Rheobase

%% Extract 1xRheobase and 2xRheobase traces 
for i=1:size(Ephys,2)
if isempty(Ephys(i).Rheobase)==0
idx_1xRheo(i)=find(Ephys(i).Rheobase.stimvec==1*Ephys(i).Rheobase.rheo);
trace_1xRheo(:,i)=Ephys(i).Rheobase.traces(:,idx_1xRheo(i));
if idx_1xRheo(i)<9
idx_2xRheo(i)=find(Ephys(i).Rheobase.stimvec==2*Ephys(i).Rheobase.rheo);
else
idx_2xRheo(i)=find(Ephys(i).Rheobase.stimvec==Ephys(i).Rheobase.stimvec(end));   
end
trace_2xRheo(:,i)=Ephys(i).Rheobase.traces(:,idx_2xRheo(i));
else
idx_1xRheo(i)=NaN;
trace_1xRheo(:,i)=ones(25000,1)*NaN;
idx_2xRheo(i)=NaN;
trace_2xRheo(:,i)=ones(25000,1)*NaN;
end
end
%% Calculate ISI between the first two spikes for 1x and 2x Rheo
for i=1:size(Ephys,2)
sp_time_2x=spike_times(trace_2xRheo(:,i),0.1);
sp_time_1x=spike_times(trace_1xRheo(:,i),0.1);
try
F2xRheo(:,i)=1/((sp_time_2x(2)-sp_time_2x(1))/20000);
catch
F2xRheo(:,i)=NaN;
end
try 
F1xRheo(:,i)=1/((sp_time_1x(2)-sp_time_1x(1))/20000);
catch
F1xRheo(:,i)=NaN;
end
end    
%% 
fig2=figure;set(fig2, 'Position', [200, -200, 1200, 1200]);set(gcf,'color','w');
 for i=1:size(Ephys,2)
hold on;subplot(7,6,i);
plot(trace_1xRheo(:,i),'Color','k','LineWidth',1);set(gca,'box','off');axis off;
hold on;plot(trace_2xRheo(:,i),'Color','r','LineWidth',1);set(gca,'box','off');axis off;
 end

%% 
dt = 1e-4;
dy = 1e-3;
props = struct('spike_finder', 2, 'threshold', 0);
traces_analysis=trace(trace_2xRheo,dt, dy, 'Analysis', props);
alltrace_info = getProfileAllSpikes(traces_analysis);
parameters_2xRheo=alltrace_info.spikes_db.data; 
F2xRheo=1/((parameters_2xRheo(2,18)-parameters_2xRheo(1,18))/2000);

%% 

                        dt = 1e-4;
                        dy = 1e-3;
                        props = struct('spike_finder', 2, 'threshold', 0);
                        traces_analysis=trace(trace_curr,dt, dy, 'Analysis', props);
                        alltrace_info = getProfileAllSpikes(traces_analysis);
                        parameters=alltrace_info.spikes_db.data; 
%% 

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


%% LED plotting


%Plotting
                       fig5= figure;set(fig5, 'Name', 'LED stimulation');set(fig5, 'Position', [200, 100, 400, 200]);set(gcf,'color','w');
                       plot(filt_traces,'Color','k','LineWidth',1);box off;axis off;ylim([min(filt_traces)-50 abs(max(filt_traces))*1.5]);
                       hold on;co=[0:1:freq-1];
                       for h=1:freq
                       x1=delay+duration/freq*co(h);
                       x2=(delay+pulsedur)+duration/freq*co(h);
                       p1=plot([x1 x2],[abs(max(filt_traces))*1.1 abs(max(filt_traces))*1.1],'-','Color','b','LineWidth',2);
                       end
                       hold on;text(0,abs(max(filt_traces))*1.5,'472 nm','Color','b');
                       scale_x= 200;
                       scale_y= 40;
                       ov_min=min(min(filt_traces));
                       %scale barx
                       hold on;x1= length(filt_traces)-200*srF;x2=length(filt_traces);p1=plot([x1 x2],[ov_min-45 ov_min-45],'-','Color','k','LineWidth',1.5);
                       %scale bary
                       hold on;y2= (ov_min-45)+scale_y;y1=(ov_min-45);p2=plot([x2 x2],[y1 y2],'-','Color','k','LineWidth',1.5); 
                       title(cell_ab);
                       
                       
                       
                       
                            %Plotting
%                        fig6= figure;set(fig6, 'Name', 'LED stimulation');set(fig6, 'Position', [200, 100, 400, 200]);set(gcf,'color','w');
%                        plot(filt_traces,'Color','k','LineWidth',1);box off;axis off;ylim([min(filt_traces)-50 abs(max(filt_traces))*1.5]);
%                        hold on;
%                        x1=delay;
%                        x2=delay+pulsedur;
%                        p1=plot([x1 x2],[abs(max(filt_traces))*1.1 abs(max(filt_traces))*1.1],'-','Color','b','LineWidth',5);
%                        hold on;text(0,abs(max(filt_traces))*1.5,'472 nm','Color','b');
%                         scale_x= 200;
%                        scale_y= abs(min(filt_traces))/10;
%                        ov_min=min(min(filt_traces));
%                        %scale barx
%                        hold on;x1= length(filt_traces)-50*srF;x2=length(filt_traces);p1=plot([x1 x2],[ov_min-45 ov_min-45],'-','Color','k','LineWidth',1.5);
%                        %scale bary
%                        hold on;y2= (ov_min-45)+scale_y;y1=(ov_min-45);p2=plot([x2 x2],[y1 y2],'-','Color','k','LineWidth',1.5); 
%                        title(cell_ab);
%                     elseif contains(stimuli_type,'pulse')==1 
%                         delay=str2num(data.header.StimulusLibrary.Stimuli.element11.Delegate.Delay)*sr;
%                        pulsedur=str2num(data.header.StimulusLibrary.Stimuli.element11.Delegate.Duration)*sr;
%                        endtime=data.header.StimulusLibrary.Stimuli.element11.Delegate.EndTime*sr;
%                        %Plotting
%                        fig7= figure;set(fig7, 'Name', 'LED stimulation');set(fig7, 'Position', [200, 100, 400, 200]);set(gcf,'color','w');
%                        plot(filt_traces,'Color','k','LineWidth',1);box off;axis off;ylim([min(filt_traces)-50 abs(max(filt_traces))*1.5]);
%                        hold on;
%                        x1=delay;
%                        x2=delay+pulsedur;
%                        p1=plot([x1 x2],[abs(max(filt_traces))*1.1 abs(max(filt_traces))*1.1],'-','Color','b','LineWidth',5);
%                        hold on;text(0,abs(max(filt_traces))*1.5,'472 nm','Color','b');
%                         scale_x= 200;
%                        scale_y= abs(min(filt_traces))/10;
%                        ov_min=min(min(filt_traces));
%                        %scale barx
%                        hold on;x1= length(filt_traces)-50*srF;x2=length(filt_traces);p1=plot([x1 x2],[ov_min-45 ov_min-45],'-','Color','k','LineWidth',1.5);
%                        %scale bary
%                        hold on;y2= (ov_min-45)+scale_y;y1=(ov_min-45);p2=plot([x2 x2],[y1 y2],'-','Color','k','LineWidth',1.5); 
%                        title(cell_ab);


 %Plotting spike frequency
%                        fig2=figure;set(fig2, 'Position', [200, 800, 200, 200]);set(gcf,'color','w');
%                        plot(stimvec,spikecount,'--bo');set(gca,'box','off');set(gcf,'color','w');xlabel('Current (pA)');ylabel('Spike Frequency (Hz)');


%Rhoebase

  
%                        %Plotting
%                        fig4=figure;set(fig4, 'Position', [200, 600, 200, 300]);set(gcf,'color','w');
%                        subplot(2,1,1);
%                        plot(traces_deconc(:,spike_idx(1)),'Color','k','LineWidth',2);set(gca,'box','off');xlim([-260*srF length(traces_deconc)]);
%                        hold on;text(-260*srF,round(parameters(3),1),[num2str(round(parameters(3),1)),'mV'],'FontSize',9);axis off;
%                        %Scale bar
%                        scale_x= 200;
%                        scale_y= 40;
%                        ov_min=min(min(traces_deconc(:,spike_idx(1))))-40;
%                        %scale barx
%                        hold on;x1= 1250*srF-(100*srF);x2=1250*srF;p1=plot([x1 x2],[ov_min-10 ov_min-10],'-','Color','k','LineWidth',1.5);
%                        %scale bary
%                        hold on;y2= (ov_min-10)+scale_y;y1=(ov_min-10);p2=plot([x2 x2],[y1 y2],'-','Color','k','LineWidth',1.5); 
%                        title(cell_ab);
%                        
%                        subplot(2,1,2);
%                        ov_min=min(min(squarepulse));
%                        plot(squarepulse(:,spike_idx(1)),'Color',[0.7 0.7 0.7],'LineWidth',1);set(gca,'box','off');xlim([-260*srF length(traces_deconc)]);
%                        axis off;
%                        %Scale bar
%                        scale_x= 200;
%                        scale_y= 20;
%                        %scale barx
%                        hold on;x1= (1250*srF)-(100*srF);x2=1250*srF;p1=plot([x1 x2],[ov_min-30 ov_min-30],'-','Color','k','LineWidth',1.5);
%                        %scale bary
%                        hold on;y2= (ov_min-30)+scale_y;y1=(ov_min-30);p2=plot([x2 x2],[y1 y2],'-','Color','k','LineWidth',1.5); 
%                        hold on;text(-260*srF,Rheobase,[num2str(round(Rheobase,1)),'pA'],'FontSize',9);axis off;

% Passive
% %Plotting fit
%                            fig3=figure;set(fig3, 'Position', [400, 800, 200, 200]);set(gcf,'color','w');
%                            plot(stimvec,dV_IR,'ko');
%                            hold on;ylabel('\DeltaVoltage (mV)');xlabel('Current (pA)');
%                            plot(stimvec,yfit,'r-.');set(gca,'box','off');axis square;

% 
% %Input resistance
%                            Rin=round(P(1)*1000,0);
%                            hold on;
%                            text(-30,max(dV_IR),[num2str(Rin),'M\Omega']);