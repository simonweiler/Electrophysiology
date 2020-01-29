
%L6 plotting scripts 
%% Load L6 structure 
str_L6    = 'C:\Users\Simon-localadmin\Documents\MargrieLab\MatlabCode\Electrophysiology\Ephys_slice\output_structure';
folder_list = uipickfiles('FilterSpec',str_L6);
load(char(folder_list));
srF=20;
sr=20000;
%% Intrinsic properties
%% Read out overall max and min + label category for figure 
step=2;
for i=1:size(Ephys,2)
    label(i)=Ephys(i).label;
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
fig1=figure;set(fig1, 'Position', [200, -200, 1400, 1000]);set(gcf,'color','w');
for i=1:size(Ephys,2)
hold on;subplot(8,9,i);
if isempty(Ephys(i).IV)==0
%Plot only every second trace atm for display
plot(Ephys(i).IV.traces(:,1:step:end),'Color','k','LineWidth',1);set(gca,'box','off');axis off;
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
% hold on;subplot(8,9,size(Ephys,2)+1);
% for m=1:length(Ephys(1).IV.stimvec);  
% squarepulse=[repmat(0,2500,1) ;repmat(Ephys(1).IV.stimvec(m),20000,1) ;repmat(0,2500,1)]
% plot(squarepulse,'Color',[0.7 0.7 0.7],'LineWidth',1);set(gca,'box','off');axis off
% hold on;
% end
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
rang=Ephys(i).Rheobase.traces(:,idx_1xRheo(i):idx_2xRheo(i));
Rheo_all{:,i}=rang;
rang_idx=[idx_1xRheo(i):idx_2xRheo(i)];
rang_all{:,i}=rang_idx;
curr_inj=Ephys(i).Rheobase.stimvec(rang_idx);
curr_all{:,i}=curr_inj;
else
idx_1xRheo(i)=NaN;
trace_1xRheo(:,i)=ones(25000,1)*NaN;
idx_2xRheo(i)=NaN;
trace_2xRheo(:,i)=ones(25000,1)*NaN;
Rheo_all{:,i}=ones(25000,1)*NaN;
rang_idx=NaN;
rang_all{:,i}=NaN;
curr_all{:,i}=NaN;
end
end
%% Calculate ISI between the first two spikes for 1x and 2x Rheo and 2xRheo 200 ms see Velez Fort 2014
delay=0.125;
spt=(delay*sr)+0.2*sr;
for i=1:size(Ephys,2)
sp_time_2x=spike_times(trace_2xRheo(:,i),0.1);
sp_time_1x=spike_times(trace_1xRheo(:,i),0.1);
spt_200=find(sp_time_2x>=spt);
try
F2xRheo(:,i)=1/((sp_time_2x(2)-sp_time_2x(1))/sr);
F2xRheo_av(:,i)=length(sp_time_2x);
catch
F2xRheo(:,i)=0;
F2xRheo_av(:,i)=0;
end
try 
F1xRheo(:,i)=1/((sp_time_1x(2)-sp_time_1x(1))/sr);
catch
F1xRheo(:,i)=0;
end
try
F2xRheo200(:,i)=1/((sp_time_2x(spt_200(2))-sp_time_2x(spt_200(1)))/sr);
catch
F2xRheo200(:,i)=0;
end
end 
%Eacc Index: early accomodation index to assess adaptation  
Eaci=((F2xRheo-F2xRheo200)./F2xRheo)*100;
Eaci(find(isnan(Eaci)))=0;
%% Extract the slope F1 Rheo1x to 2xRheo vs injected current
for i=1:size(Ephys,2)
for k=1:size(Rheo_all{:,i},2)
sp_time_all=spike_times(Rheo_all{:,i}(:,k),0.1);
try
F_idm{i,k}=1/((sp_time_all(2)-sp_time_all(1))/sr);
F_h=1/((sp_time_all(2)-sp_time_all(1))/sr);
catch
 F_idm{i,k}=NaN;   
end
end
end
%% Slope for injected current vs inst.fre from 1xRheo to 2xRheo
for i=1:size(Ephys,2)
  ifreq=[F_idm{i,:}];
  curri=curr_all{1,i};
  P = polyfit(curri(find(~isnan(ifreq)))', ifreq(find(~isnan(ifreq)))',1);
  F1_slope(i)=P(1); 
end
%% Plot the 1x and 2x Rheobase traces
ov_min=nanmin(nanmin([trace_1xRheo trace_2xRheo]));
ov_max=nanmax(nanmax([trace_1xRheo trace_2xRheo]));
fig2=figure;set(fig2, 'Position', [200, -200, 1400, 1000]);set(gcf,'color','w');
 for i=1:size(Ephys,2)
hold on;subplot(8,9,i);
plot(trace_1xRheo(:,i),'Color','k','LineWidth',1);set(gca,'box','off');axis off;
hold on;plot(trace_2xRheo(:,i),'Color','r','LineWidth',1);set(gca,'box','off');axis off;
hold on;ylim([ov_min-10 ov_max]);
%title([num2str(F1xRheo(i)) 'Hz , ' num2str(F2xRheo(i)) 'Hz , ' num2str(Ephys(i).label)]);
title([num2str(F1xRheo(i)) ', ' num2str(F2xRheo(i)) 'Hz']);
 end
% legend('1xRheo', '2xRheo');
% legend boxoff;
%% Plot the 3 parameters: F1/I slope, Eacc index and F2Rheo against each other
fig5= figure;set(fig5, 'Name', 'Plot3');set(fig5, 'Position', [200, 100, 800, 400]);set(gcf,'color','w');
subplot(1,2,1);
for i=1:size(Ephys,2)
    if isempty(Ephys(i).IV)==0
    plot(Ephys(i).IV.stimvec,Ephys(i).IV.spikecount,'--ok'); 
    hold on;
    else       
    end
end
xlim([0 Ephys(1).IV.stimvec(end)]);
ylabel('Spike frequency (Hz');xlabel('Current (pA)');box off
subplot(1,2,2);plot3(F2xRheo,F1_slope,Eaci,'o','Color','k','MarkerFaceColor','k');grid on;set(gca,'xdir','reverse');
xlabel('F12xRb');ylabel('F1/I slope');zlabel('Eacc index');
%% Plot the 3 parameters: F1/I slope, Eacc index and F2Rheo against each other + color coded for labelled cells
pointsize=20;
fig5= figure;set(fig5, 'Name', 'Plot3');set(fig5, 'Position', [200, 100, 400, 400]);set(gcf,'color','w');
scatter3(F2xRheo,F1_slope,Eaci,pointsize,label,'filled');grid on;set(gca,'xdir','reverse');
xlabel('F12xRb');ylabel('F1/I slope');zlabel('Eacc index');
 %% Plot the ramp traces
 fig3=figure;set(fig3, 'Position', [200, -200, 1400, 1000]);set(gcf,'color','w');
for i=1:size(Ephys,2)
hold on;subplot(8,9,i);
if isempty(Ephys(i).Ramp)==0
%Plot only every second trace atm for display
plot(Ephys(i).Ramp.traces,'Color','k','LineWidth',1);set(gca,'box','off');axis off;
%xlim([-260*srF length(traces_deconc)]);
% hold on;ylim([ov_min-10 ov_max]);
% hold on;text(-550*srF,Ephys(i).IV.RMP,[num2str(Ephys(i).IV.RMP),'mV'],'FontSize',9);
 %Scale bar
%  scale_x= 200;
%  scale_y= 40;
%  %scale barx
%  hold on;x1= (1250*srF)-(100*srF);x2=1250*srF;p1=plot([x1 x2],[ov_min-10 ov_min-10],'-','Color','k','LineWidth',1.5);
%  %scale bary
%  hold on;y2= (ov_min-10)+scale_y;y1=(ov_min-10);p2=plot([x2 x2],[y1 y2],'-','Color','k','LineWidth',1.5); 
title([num2str(Ephys(i).patching_date) Ephys(i).celID]);
else
plot(1,1);set(gca,'box','off');axis off;
end
end
%%  %% Plot the sag traces
 fig10=figure;set(fig10, 'Position', [200, -200, 1400, 1000]);set(gcf,'color','w');
for i=1:size(Ephys,2)
hold on;subplot(8,9,i);
if isempty(Ephys(i).Sag)==0
%Plot only every second trace atm for display
plot(Ephys(i).Sag.traces(:,:),'Color','k','LineWidth',1);set(gca,'box','off');axis off;
%xlim([-260*srF length(traces_deconc)]);
% hold on;ylim([ov_min-10 ov_max]);
% hold on;text(-550*srF,Ephys(i).IV.RMP,[num2str(Ephys(i).IV.RMP),'mV'],'FontSize',9);
 %Scale bar
%  scale_x= 200;
%  scale_y= 40;
%  %scale barx
%  hold on;x1= (1250*srF)-(100*srF);x2=1250*srF;p1=plot([x1 x2],[ov_min-10 ov_min-10],'-','Color','k','LineWidth',1.5);
%  %scale bary
%  hold on;y2= (ov_min-10)+scale_y;y1=(ov_min-10);p2=plot([x2 x2],[y1 y2],'-','Color','k','LineWidth',1.5); 
title([num2str(Ephys(i).patching_date) Ephys(i).celID]);
else
plot(1,1);set(gca,'box','off');axis off;
end
end

%% Gather parameters
%Parameters are:
% %%%%[ 1]    'MinVm'        
%     [ 2]    'PeakVm'       
%     [ 3]    'InitVm'       
%     [ 4]    'InitVmBySlope'
%     [ 5]    'MaxVmSlope'   
%     [ 6]    'HalfVm'       
%     [ 7]    'Amplitude'    
%     [ 8]    'MaxAHP'       
%     [ 9]    'DAHPMag'      
%     [10]    'InitTime'     
%     [11]    'RiseTime'     
%     [12]    'FallTime'     
%     [13]    'MinTime'      
%     [14]    'BaseWidth'    
%     [15]    'HalfWidth'    
%     [16]    'FixVWidth'    
%     [17]    'Index'        
%     [18]    'Time'
for i=1:size(Ephys,2)
    if isempty(Ephys(i).Rheobase)==0
   Vmin(i) = Ephys(i).Rheobase.parameters(1,1);
   Vpeak(i) = Ephys(i).Rheobase.parameters(1,2);
   Vthresh(i) = Ephys(i).Rheobase.parameters(1,3);
   Vslope(i) = Ephys(i).Rheobase.parameters(1,5);
   Vhalf(i) = Ephys(i).Rheobase.parameters(1,6);
   Vamp(i) = Ephys(i).Rheobase.parameters(1,7);
   Vmaxahp(i)=Ephys(i).Rheobase.parameters(1,8);
   Vinitime(i)=Ephys(i).Rheobase.parameters(1,10);
   Vrisetime(i)=Ephys(i).Rheobase.parameters(1,11);
   Vfalltime(i)=Ephys(i).Rheobase.parameters(1,12);
   Vhwidth(i)=Ephys(i).Rheobase.parameters(1,15);
   freq_max(i)=max(Ephys(i).IV.spikecount);
   rin(i)=Ephys(i).Passive.Rin;
   tau(i)=nanmean(Ephys(i).Passive.tau);
%    if Ephys(i).sag==1
%       sag(i)=Ephys(i).Sag.sagratio(1);
      
   %Ephys(i).sag==0
       sagtrace=Ephys(i).IV.traces(:,1);
       Vrest_sag=mean(sagtrace(1:delay*sr,:));
       Vss_sag=mean(sagtrace(7*delay*sr:8.9*delay*sr,:));
       Vmin_sag=min(sagtrace(delay*sr:2*delay*sr,:));
      sag(i)=100*((Vss_sag-Vmin_sag)./(Vrest_sag-Vmin_sag));
          
    else
    Vmin(i) = NaN;
   Vpeak(i) = NaN;
   Vthresh(i) = NaN;
   Vslope(i) = NaN;
   Vhalf(i) = NaN;
   Vamp(i) = NaN;
   Vmaxahp(i)=NaN;
   Vinitime(i)=NaN;
   Vrisetime(i)=NaN;
   Vfalltime(i)=NaN;
   Vhwidth(i)=NaN;
   freq_max(i)=NaN;
   rin(i)=NaN;
   tau(i)=NaN;
    sag(i)=NaN;
    end
   
end
%% 
 %param=[F2xRheo' F2xRheo_av' F1_slope' Eaci' freq_max' rin' tau' sag' Vmin' Vpeak' Vthresh' Vslope' Vhalf' Vamp'...
    % Vmaxahp' Vinitime' Vrisetime' Vfalltime' Vhwidth' label'];

 param=[rin' tau' sag' Vmin' Vpeak' Vthresh' Vslope' Vhalf' Vamp'...
     Vmaxahp' Vinitime' Vrisetime' Vfalltime' Vhwidth' ];

%% Plot histograms of all intrinsic properties
fig7=figure;set(fig7, 'Position', [200, 200, 1400, 1000]);set(gcf,'color','w');
%idx_u=find(sum(~isnan(param),2)==size(param,2));
for i=1:size(param,2)
    subplot(5,4,i)
    histogram(param(find(label==0),i),10,'FaceColor','k');box off;ylabel('Counts');
    hold on;
    histogram(param(find(label==1),i),10,'FaceColor','b');box off;ylabel('Counts');
end
%% Plot histogram of selected parameters
fig8=figure;set(fig8, 'Position', [200, 200, 200, 200]);set(gcf,'color','w');
histogram(Eaci(find(label<2)),5,'FaceColor','k');box off;
xlabel('Eacc index');ylabel('Counts');
%% Plot the parametrs
% str={'V_{min}(mV)','V_{peak}(mV)','V_{thresh}(mV)', 'Vslope_{max} (\DeltamV/\Deltams)','V_{half} (mV)','Spike_{amplitude} (mV)',...
%     'AHP_{max}(mV)','Spike_{rise} (ms)', 'Spike_{fall} (ms)','Spike_{base width}(ms)', 'Spike_{half width} (ms)','First spike latency (ms)',...
%     'V_{rest} (mV)','Tau (ms)','R_{IN}(M\Omega)','Sag ratio (%)','Rheobase (pA)','Spike frequency_{max} (Hz)','Pial depth (µm)'};
str={'R_{IN}(M\Omega)','Tau (ms)','Sag ratio (%)','V_{min}(mV)','V_{peak}(mV)','V_{thresh}(mV)', 'Vslope_{max} (\DeltamV/\Deltams)','V_{half} (mV)','Spike_{amplitude} (mV)',...
    'AHP_{max}(mV)', 'Spike init (ms)','Spike_{rise} (ms)', 'Spike_{fall} (ms)','Spike_{half width} (ms)'};

fig9=figure;set(fig9, 'Position', [200, 200, 1400, 400]);set(gcf,'color','w');
idx_u=find(sum(~isnan(param),2)==size(param,2));
for i=1:size(param,2)
    subplot(2,7,i)
    histogram(param(idx_u,i),10,'FaceColor','k');box off;ylabel('Counts');
    hold on;
    title(str{i});
end


%% Clustering
clu_num = 3;
pia=rand(length(param(idx_u,:)),1);
[idx_input, clustering_input, leafOrder] = hca(param(idx_u,:),1,'ward',clu_num,pia,1);%call function for clustering

%% Plot train postive current clamp
 fig3=figure;set(fig3, 'Position', [200, -200, 1400, 1000]);set(gcf,'color','w');
 for i=1:size(Ephys,2)
hold on;subplot(7,7,i);
if size(Ephys(i).sub_traces_train,2)>1
  plot(Ephys(i).sub_traces_train(:,1),'Color','k','LineWidth',1);set(gca,'box','off');%axis off;  
else
 plot(1,1);set(gca,'box','off');axis off;   
end
 end
 %% 

% dt = 1e-4;
% dy = 1e-3;
% props = struct('spike_finder', 2, 'threshold', 0);
% traces_analysis=trace(trace_2xRheo,dt, dy, 'Analysis', props);
% alltrace_info = getProfileAllSpikes(traces_analysis);
% parameters_2xRheo=alltrace_info.spikes_db.data; 
% F2xRheo=1/((parameters_2xRheo(2,18)-parameters_2xRheo(1,18))/2000);

%% 

%                         dt = 1e-4;
%                         dy = 1e-3;
%                         props = struct('spike_finder', 2, 'threshold', 0);
%                         traces_analysis=trace(trace_curr,dt, dy, 'Analysis', props);
%                         alltrace_info = getProfileAllSpikes(traces_analysis);
%                         parameters=alltrace_info.spikes_db.data; 
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