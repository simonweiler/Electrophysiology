
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
    if isempty(Ephys(i).IV)==0;
    rmp(i)=Ephys(i).IV.RMP;
    else
    rmp(i)=NaN;
    end
    %pia(i)=Ephys(i).pia;
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
step=3
fig1=figure;set(fig1, 'Position', [200, -200, 1400, 1000]);set(gcf,'color','w');
for i=1:size(Ephys,2)
hold on;subplot(10,11,i);
if isempty(Ephys(i).IV)==0
%Plot only every second trace atm for display
if label(i)==1
plot(Ephys(i).IV.traces(:,1:step:end),'Color','g','LineWidth',1);set(gca,'box','off');axis off;
elseif label(i)==0
plot(Ephys(i).IV.traces(:,1:step:end),'Color','k','LineWidth',1);set(gca,'box','off');axis off;
else
plot(Ephys(i).IV.traces(:,1:step:end),'Color','m','LineWidth',1);set(gca,'box','off');axis off;
end
%xlim([-260*srF length(traces_deconc)]);
hold on;ylim([ov_min-10 ov_max]);
hold on;text(-550*srF,Ephys(i).IV.RMP,[num2str(rmp(i)),'mV'],'FontSize',9);
 %Scale bar
 scale_x= 200;
 scale_y= 40;
 %scale barx
 hold on;x1= (1250*srF)-(100*srF);x2=1250*srF;p1=plot([x1 x2],[ov_min-10 ov_min-10],'-','Color','k','LineWidth',1.5);
 %scale bary
 hold on;y2= (ov_min-10)+scale_y;y1=(ov_min-10);p2=plot([x2 x2],[y1 y2],'-','Color','k','LineWidth',1.5); 
title([num2str(Ephys(i).patching_date) Ephys(i).celID ]);
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
%% Read out labeled and non labelled/ INs
la=find(label==1);
nla=find(label==0);
in=find(label==2);
%% Plot Frequcny 
fig5= figure;set(fig5, 'Name', 'Plot3');set(fig5, 'Position', [200, 100, 800, 400]);set(gcf,'color','w');
subplot(1,2,1);
for i=1:size(Ephys,2)
    if isempty(Ephys(i).IV)==0
    plot(Ephys(i).IV.stimvec,Ephys(i).IV.spikecount,'--ok'); 
    max_spikefre(i)=max(Ephys(i).IV.spikecount);
    hold on;
    else  
      max_spikefre(i)=NaN;  
    end
end
xlim([0 Ephys(1).IV.stimvec(end)]);
ylabel('Spike frequency (Hz');xlabel('Current (pA)');box off
%% Hisogram of max spikefrequency all 3 groups
figure;set(gcf,'color','w');histogram(max_spikefre(la),10,'FaceColor','g');hold on;histogram(max_spikefre(nla),10,'FaceColor','k');hold on;histogram(max_spikefre(in),15,'FaceColor','m');box off;
legend('CPN','n-CPN','n-CPN IN');xlabel('max AP frequency (Hz)');ylabel('Cells');legend boxoff ;
%% Hisogram of max spikefrequency all 2 groups
figure;set(gcf,'color','w');histogram(max_spikefre(la),14,'FaceColor','g');hold on;histogram(max_spikefre(nla),14,'FaceColor','k');box off;
legend('CPN','n-CPN');xlabel('max AP frequency (Hz)');ylabel('Cells');legend boxoff ;
%% %% Gather parameters
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
%    freq_max(i)=max(Ephys(i).IV.spikecount);
%    rin(i)=Ephys(i).Passive.Rin;
%    tau(i)=nanmean(Ephys(i).Passive.tau);
%    if Ephys(i).sag==1
%       sag(i)=Ephys(i).Sag.sagratio(1);
      
   %Ephys(i).sag==0
      
          
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
  
    %sag(i)=NaN;
    end
   
end
%% 
for i=1:size(Ephys,2)
 if isempty(Ephys(i).Passive)==0
 rin(i)=Ephys(i).Passive.Rin;
 tau(i)=nanmean(Ephys(i).Passive.tau)
 else
   rin(i)=NaN;
   tau(i)=NaN;   
    end
    end
%% 
for i=1:size(Ephys,2)
 if isempty(Ephys(i).Sag)==0
     if Ephys(i).sag==1
         if Ephys(i).Sag.sagratio(1)>0
            sag(i)=Ephys(i).Sag.sagratio(1);
         else
            sag(i)=Ephys(i).Sag.sagratio(2);
         end
     else
sag(i)=NaN;
     end
 else
 sag(i)=NaN;
 end
end
%  Vrest_sag=mean(sagtrace(1:delay*sr,:));
%        Vss_sag=mean(sagtrace(7*delay*sr:8.9*delay*sr,:));
%        Vmin_sag=min(sagtrace(delay*sr:2*delay*sr,:));
%       sag(i)=100*((Vss_sag-Vmin_sag)./(Vrest_sag-Vmin_sag));
%% 
 %param=[F2xRheo' F2xRheo_av' F1_slope' Eaci' freq_max' rin' tau' sag' Vmin' Vpeak' Vthresh' Vslope' Vhalf' Vamp'...
    % Vmaxahp' Vinitime' Vrisetime' Vfalltime' Vhwidth' label'];

 param=[sag' Vmin' Vpeak' Vthresh' Vslope' Vhalf' Vamp'...
     Vmaxahp' (Vinitime/2)' (Vrisetime/2)' (Vfalltime/2)' (Vhwidth/2)' ];

%% Plot histograms of all intrinsic properties
fig9=figure;set(fig9, 'Position', [200, 200, 1400, 400]);set(gcf,'color','w');
%idx_u=find(sum(~isnan(param),2)==size(param,2));
str={'Sag ratio (%)','V_{min}(mV)','V_{peak}(mV)','V_{thresh}(mV)', 'Vslope_{max} (\DeltamV/\Deltams)','V_{half} (mV)','Spike_{amplitude} (mV)',...
    'AHP_{max}(mV)', 'Spike init (ms)','Spike_{rise} (ms)', 'Spike_{fall} (ms)','Spike_{half width} (ms)'};

for i=1:size(param,2)
    subplot(2,6,i)
     histogram(param(find(label==1),i),10,'FaceColor','g');box off;ylabel('Cells');
     hold on;
    histogram(param(find(label==0),i),10,'FaceColor','k');box off;ylabel('Cells');
    hold on;
   
   % histogram(param(find(label==2),i),10,'FaceColor','m');box off;ylabel('Cells');
     hold on;
    title(str{i});
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


%% Rheobase
trace_1xRheo_2ap=[];
for i=1:size(Ephys,2)
    if isempty(Ephys(i).Rheobase)==0;
        %find Rheobase 
        rheo(i)=Ephys(i).Rheobase.rheo;
        idx_1xRheo(i)=find(Ephys(i).Rheobase.stimvec==1*Ephys(i).Rheobase.rheo);
        trace_1xRheo(:,i)=Ephys(i).Rheobase.traces(:,idx_1xRheo(i));
        trace_1xRheo_2ap(:,i)=Ephys(i).Rheobase.traces(:,idx_1xRheo(i)+2);
        sp_time_1x=spike_times(trace_1xRheo(:,i),0.1);
        %idx_xxRheo(i)=
    else
        rheo(i)=NaN;
        idx_1xRheo(i)=NaN;
        trace_1xRheo(:,i)=ones(length(Ephys(1).Rheobase.traces),1)*NaN;
    end
end
%% 
fig1=figure;set(fig1, 'Position', [200, -200, 1400, 1000]);set(gcf,'color','w');
hold on;
for i=1:size(Ephys,2)
    subplot(10,11,i);
   plot(trace_1xRheo_2ap(:,i)); box off;
end

%% Plot  of Rheo and rmp
par=rin;
yl='R_{IN}(M\Omega)';
%yl='Rheobase (pA)';
comi=[par(la) par(nla) par(in)];
group=[ones(1,length(par(la))) ones(1,length(par(nla)))*2 ones(1,length(par(in)))*3];
figure; set(gcf,'color','w');
plot(group(1:length(par(la))),par(la),'go');hold on;plot(ones(1,length(par(nla)))*2,par(nla),'ko');hold on;plot(ones(1,length(par(in)))*3,par(in),'mo')
%plotSpread(comi'-90,'categoryIdx',group, 'categoryMarkers',{'.','.'},'categoryColors',{'k','m'});
hold on;
boxplot(comi,group,'Colors','k','OutlierSize',0.001);
box off;
xticklabels({'CPNs','n-CPNs','IN'});ylabel(yl);
%% Plot three examples for CPN, n-CPN and IN
fig4=figure;set(fig4, 'Position', [200, 800, 600, 200]);set(gcf,'color','w');
subplot(1,3,1)
plot(Ephys(30).IV.traces(:,1:step:end),'Color','g','LineWidth',1);set(gca,'box','off');axis off;
hold on;ylim([ov_min-10 ov_max]);
hold on;text(-550*srF,Ephys(30).IV.RMP,[num2str(rmp(30)),'mV'],'FontSize',9);
hold on;title('CPN');
subplot(1,3,2)
plot(Ephys(67).IV.traces(:,1:step:end),'Color','k','LineWidth',1);set(gca,'box','off');axis off;
hold on;ylim([ov_min-10 ov_max]);
hold on;text(-550*srF,Ephys(67).IV.RMP,[num2str(rmp(67)),'mV'],'FontSize',9);hold on;title('n-CPN');
subplot(1,3,3)
plot(Ephys(15).IV.traces(:,1:step:end),'Color','m','LineWidth',1);set(gca,'box','off');axis off;
hold on;ylim([ov_min-10 ov_max]);
hold on;text(-550*srF,Ephys(15).IV.RMP,[num2str(rmp(15)),'mV'],'FontSize',9);
hold on;title('n-CPN FS-IN');

%% Extract 1xRheobase and 2xRheobase traces 
trace_xxRheo=[];
idx_xxRheo=[];
for i=1:size(Ephys,2)
if isempty(Ephys(i).Rheobase)==0
idx_1xRheo(i)=find(Ephys(i).Rheobase.stimvec==1*Ephys(i).Rheobase.rheo);
trace_1xRheo(:,i)=Ephys(i).Rheobase.traces(:,idx_1xRheo(i));
if idx_1xRheo(i)<9
idx_2xRheo(i)=find(Ephys(i).Rheobase.stimvec==2*Ephys(i).Rheobase.rheo);
else
idx_2xRheo(i)=find(Ephys(i).Rheobase.stimvec==Ephys(i).Rheobase.stimvec(end));  
 
end
if idx_1xRheo(i)<7
idx_xxRheo(i)=find(Ephys(i).Rheobase.stimvec==round(2.7*Ephys(i).Rheobase.rheo,-1));
else
idx_xxRheo(i)=find(Ephys(i).Rheobase.stimvec==Ephys(i).Rheobase.stimvec(end)); 
end
%extract end-3
%idx_xxRheo(i)=find(Ephys(i).Rheobase.stimvec==Ephys(i).Rheobase.stimvec(end-1));  
trace_2xRheo(:,i)=Ephys(i).Rheobase.traces(:,idx_2xRheo(i));
trace_xxRheo(:,i)=Ephys(i).Rheobase.traces(:,idx_xxRheo(i));
rang=Ephys(i).Rheobase.traces(:,idx_1xRheo(i):idx_2xRheo(i));
rangxx=Ephys(i).Rheobase.traces(:,idx_1xRheo(i):idx_xxRheo(i));
Rheo_all{:,i}=rang;
Rheo_allxx{:,i}=rangxx;
rang_idx=[idx_1xRheo(i):idx_2xRheo(i)];
rang_idx_xx=[idx_1xRheo(i):idx_xxRheo(i)];
rang_all{:,i}=rang_idx;
rang_allxx{:,i}=rang_idx_xx;
curr_inj=Ephys(i).Rheobase.stimvec(rang_idx);
curr_injxx=Ephys(i).Rheobase.stimvec(rang_idx_xx);
curr_all{:,i}=curr_inj;
curr_allxx{:,i}=curr_injxx;
else
idx_1xRheo(i)=NaN;
trace_1xRheo(:,i)=ones(25000,1)*NaN;
idx_2xRheo(i)=NaN;
idx_xxRheo(i)=NaN;
trace_2xRheo(:,i)=ones(25000,1)*NaN;
trace_xxRheo(:,i)=ones(25000,1)*NaN;
Rheo_all{:,i}=ones(25000,1)*NaN;
Rheo_allxx{:,i}=ones(25000,1)*NaN;
rang_idx_xx=NaN;
rang_allxx{:,i}=NaN;
curr_allxx{:,i}=NaN;
end
end
%% 

plot_Rheobase(Ephys,la);
%% 

for i=1:size(Ephys,2)
if isempty(Ephys(i).Rheobase)==0
idx_1xRheo(i)=find(Ephys(i).Rheobase.stimvec==1*Ephys(i).Rheobase.rheo);
trace_1xRheo(:,i)=Ephys(i).Rheobase.traces(:,idx_1xRheo(i));
rang_all=Ephys(i).Rheobase.traces(:,idx_1xRheo(i):end);
Rheo_all_full{:,i}=rang_all;
else
idx_1xRheo(i)=NaN;
trace_1xRheo(:,i)=ones(25000,1)*NaN;
rang_all=NaN;
Rheo_all_full{:,i}=ones(25000,1)*NaN;
end
end
%% 

for i=1:size(Ephys,2)
    temp=Rheo_all_full{:,i};
    for k=1:size(temp,2)
    st=spike_times(temp(:,k),0.1);
    if length(st)>=2
    diff_st(k)=st(2)-st(1);
    elseif length(st)==1
    diff_st(k)=NaN;
    else
    diff_st(k)=NaN;    
    end
    end
    diff_all{:,i}=diff_st;
   try
    diff_shtemp=find(diff_st==nanmin(diff_st));
   diff_sh(:,i)=diff_shtemp(1);
   diff_shortest(i)=diff_st(diff_sh(i));
   catch
      diff_sh(:,i)=NaN; 
      diff_shortest(i)=NaN;
   end
 
    diff_st=[];
end
%% 

1./diff_shortest*sr
%% 


figure;hold on;
for i=1:length(la)
plot(diff_all{:,la(i)},'--go')
end
hold on;
for i=1:length(nla)
plot(diff_all{:,nla(i)},'--ro')
end

%% Calculate ISI between the first two spikes for 1x and 2x Rheo and 2xRheo 200 ms see Velez Fort 2014
delay=0.125;
spt=(delay*sr)+0.2*sr;
for i=1:size(Ephys,2)
sp_time_2x=spike_times(trace_2xRheo(:,i),0.1);
sp_time_xx=spike_times(trace_xxRheo(:,i),0.1);
sp_time_1x=spike_times(trace_1xRheo(:,i),0.1);
spt_200=find(sp_time_2x>=spt);
spt_200xx=find(sp_time_xx>=spt);
try
F2xRheo(:,i)=1/((sp_time_2x(2)-sp_time_2x(1))/sr);
F2xRheoxx(:,i)=1/((sp_time_xx(2)-sp_time_xx(1))/sr);
F2xRheo_av(:,i)=length(sp_time_2x);
F2xRheo_avxx(:,i)=length(sp_time_xx);
catch
F2xRheo(:,i)=0;
F2xRheoxx(:,i)=0;
F2xRheo_av(:,i)=0;
F2xRheo_avxx(:,i)=0;
end
try 
F1xRheo(:,i)=1/((sp_time_1x(2)-sp_time_1x(1))/sr);
catch
F1xRheo(:,i)=0;
end
try
F2xRheo200(:,i)=1/((sp_time_2x(spt_200(2))-sp_time_2x(spt_200(1)))/sr);
F2xRheo200xx(:,i)=1/((sp_time_xx(spt_200xx(2))-sp_time_xx(spt_200xx(1)))/sr);
catch
F2xRheo200(:,i)=0;
F2xRheo200xx(:,i)=0;
end
end 
%Eacc Index: early accomodation index to assess adaptation  
Eaci=((F2xRheo-F2xRheo200)./F2xRheo)*100;
Eacixx=((F2xRheoxx-F2xRheo200xx)./F2xRheoxx)*100;
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
%% Extract the slope F1 Rheo1x to xxRheo vs injected current
for i=1:size(Ephys,2)
for k=1:size(Rheo_allxx{:,i},2)
sp_time_allxx=spike_times(Rheo_allxx{:,i}(:,k),0.1);
try
F_idmxx{i,k}=1/((sp_time_allxx(2)-sp_time_allxx(1))/sr);
F_hxx=1/((sp_time_allxx(2)-sp_time_allxx(1))/sr);
catch
 F_idmxx{i,k}=NaN;   
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
%% Slope for injected current vs inst.fre from 1xRheo to xxRheo
for i=1:size(Ephys,2)
  ifreq=[F_idmxx{i,:}];
  curri=curr_allxx{1,i};
  P = polyfit(curri(find(~isnan(ifreq)))', ifreq(find(~isnan(ifreq)))',1);
  F1_slopexx(i)=P(1); 
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
%% %% Plot the 1x and xx Rheobase traces
ov_min=nanmin(nanmin([trace_1xRheo trace_xxRheo]));
ov_max=nanmax(nanmax([trace_1xRheo trace_xxRheo]));
fig2=figure;set(fig2, 'Position', [200, -200, 1400, 1000]);set(gcf,'color','w');
 for i=1:size(Ephys,2)
hold on;subplot(9,9,i);
plot(trace_1xRheo(:,i),'Color','k','LineWidth',1);set(gca,'box','off');axis off;
hold on;plot(trace_xxRheo(:,i),'Color','r','LineWidth',1);set(gca,'box','off');axis off;
hold on;ylim([ov_min-10 ov_max]);
title([num2str(F1xRheo(i)) 'Hz , ' num2str(F2xRheoxx(i)) 'Hz , ' num2str(Ephys(i).label)]);
%title([num2str(F1xRheo(i)) ', ' num2str(F2xRheoxx(i)) 'Hz']);
 end
% legend('1xRheo', '2xRheo');
% legend boxoff;

%% Plot the 3 parameters: F1/I slope, Eacc index and F2Rheo against each other

%% 

% subplot(1,2,2);plot3(F2xRheo,F1_slope,Eaci,'o','Color','k','MarkerFaceColor','k');grid on;set(gca,'xdir','reverse');
% xlabel('F12xRb');ylabel('F1/I slope');zlabel('Eacc index');
%% Plot the 3 parameters: F1/I slope, Eacc index and F2Rheo against each other + color coded for labelled cells
pointsize=20;
fig5= figure;set(fig5, 'Name', 'Plot3');set(fig5, 'Position', [200, 100, 400, 400]);set(gcf,'color','w');
scatter3(F2xRheo,F1_slope,Eaci,pointsize,label,'filled');grid on;set(gca,'xdir','reverse');
xlabel('F12xRb');ylabel('F1/I slope');zlabel('Eacc index');
%% Plot the 3 parameters: F1/I slopexx, Eaccxx index and F2Rheoxx against each other + color coded for labelled cells
pointsize=20;
fig5= figure;set(fig5, 'Name', 'Plot3');set(fig5, 'Position', [200, 100, 400, 400]);set(gcf,'color','w');
scatter3(F2xRheoxx,F1_slopexx,Eacixx,pointsize,label,'filled');grid on;set(gca,'xdir','reverse');
xlabel('F12xRb');ylabel('F1/I slope');zlabel('Eacc index');
 %% Plot the ramp traces
 fig3=figure;set(fig3, 'Position', [200, -200, 1400, 1000]);set(gcf,'color','w');
for i=1:size(Ephys,2)
hold on;subplot(9,9,i);
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



%% Clustering
clu_num = 3;
pia=rand(length(param(idx_u,:)),1);
[idx_input, clustering_input, leafOrder] = hca(param(idx_u,:),1,'ward',clu_num,pia,1);%call function for clustering










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CRACM%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot train postive current clamp
 fig3=figure;set(fig3, 'Position', [200, -200, 1400, 1000]);set(gcf,'color','w');
 for i=1:size(Ephys,2)
hold on;subplot(10,11,i);
if size(Ephys(i).sub_traces_train,2)>1
  plot(Ephys(i).sub_traces_train(:,2),'Color','k','LineWidth',1);set(gca,'box','off');%axis off;  
else
 plot(1,1);set(gca,'box','off');axis off;   
end
 end
 
 %% Read out nonzeros mV
 
 for i=1:length(la)
     if Ephys(la(i)).wc==1 & Ephys(la(i)).sol==1
     if size(Ephys(la(i)).train_p,2)>=2
    la_cc(i)=max(Ephys(la(i)).train_p(:,2));
     else
     la_cc(i)=NaN;
     end
     else
     la_cc(i)=NaN;
     end
 end
 
 
  for i=1:length(nla)
      if Ephys(nla(i)).wc==1 & Ephys(nla(i)).sol==1
     if size(Ephys(nla(i)).train_p,2)>=2
    nla_cc(i)=max(Ephys(nla(i)).train_p(:,2));
     else
     nla_cc(i)=NaN;
     end
     else
       nla_cc=NaN;
      end
  end
  
  for i=1:length(in)
      if Ephys(in(i)).wc==1 & Ephys(in(i)).sol==1
     if size(Ephys(in(i)).train_p,2)>=2
    in_cc(i)=max(Ephys(in(i)).train_p(:,2));
     else
     in_cc(i)=NaN;
     end
      else
        in_cc(i)=NaN;  
     end
  end
 %% Read  out nonzeros pA
 for i=1:length(la)
     if Ephys(la(i)).wc==1 & Ephys(la(i)).sol==1
     if size(Ephys(la(i)).train_n,2)>=1
    la_vc(i)=min(Ephys(la(i)).train_n(:,1));
     else
     la_vc(i)=NaN;
     end
     else
         la_vc(i)=NaN;
     end
     end
  
 
  for i=1:length(nla)
      if Ephys(nla(i)).wc==1 & Ephys(nla(i)).sol==1
     if size(Ephys(nla(i)).train_n,2)>=1
    nla_vc(i)=min(Ephys(nla(i)).train_n(:,1));
     else
     nla_vc(i)=NaN;
     end
      else
          nla_vc(i)=NaN;
      end
     end
 
  
     
  for i=1:length(in)
      if Ephys(in(i)).wc==1 & Ephys(in(i)).sol==1
     if size(Ephys(in(i)).train_n,2)>=1
    in_vc(i)=min(Ephys(in(i)).train_n(:,1));
     else
     in_vc(i)=NaN;
     end
      else
      in_vc(i)=NaN;
      end
  end
  %% 
   %% Read  out nonzeros pA
 for i=1:length(la)
     if Ephys(la(i)).wc==1 & Ephys(la(i)).sol==1
     if size(Ephys(la(i)).train_n,2)>=1
    la_vc_d(i)=nanmean(nonzeros(Ephys(la(i)).dpeak_n(:,1)));
     else
     la_vc_d(i)=NaN;
     end
     else
         la_vc_d(i)=NaN;
     end
     end
  

 
   for i=1:length(nla)
     if Ephys(nla(i)).wc==1 & Ephys(nla(i)).sol==1
     if size(Ephys(nla(i)).train_n,2)>=1
    nla_vc_d(i)=nanmean(nonzeros(Ephys(nla(i)).dpeak_n(:,1)));
     else
     nla_vc_d(i)=NaN;
     end
     else
         nla_vc_d(i)=NaN;
     end
     end
  
     
for i=1:length(in)
     if Ephys(in(i)).wc==1 & Ephys(in(i)).sol==1
     if size(Ephys(in(i)).train_n,2)>=1
    in_vc_d(i)=nanmean(nonzeros(Ephys(in(i)).dpeak_n(:,1)));
     else
     in_vc_d(i)=NaN;
     end
     else
         in_vc_d(i)=NaN;
     end
     end
  
 
 %% Represantative example traces for nCPN CP and IN
fig5=figure;set(fig5, 'Position', [200, 800, 600, 200]);set(gcf,'color','w');
subplot(1,3,1)
plot(Ephys(34).sub_traces_train(1:100000,2),'Color','g','LineWidth',1);set(gca,'box','off');axis off;
hold on;ylim([-1 5]);
freq=1;
co=[0:1:5-1];
duration=20000;
pulsedur=1000;
delay=5000;
  for h=1:5
  x1=delay+duration/freq*co(h);
  x2=(delay+pulsedur)+duration/freq*co(h);
  p1=plot([x1 x2],[abs(max(Ephys(34).sub_traces_train(:,2)))*1.1 abs(max(Ephys(34).sub_traces_train(:,2)))*1.1],'-','Color','b','LineWidth',2);
  hold on;
  end
%hold on;text(-550*srF,Ephys(7).IV.RMP,[num2str(rmp(7)),'mV'],'FontSize',9);
hold on;title('CPN');
subplot(1,3,2)
plot(Ephys(47).sub_traces_train(1:100000,2),'Color','k','LineWidth',1);set(gca,'box','off');axis off;
hold on;ylim([-1 5]);
  for h=1:5
  x1=delay+duration/freq*co(h);
  x2=(delay+pulsedur)+duration/freq*co(h);
  p1=plot([x1 x2],[abs(max(Ephys(34).sub_traces_train(:,2)))*1.1 abs(max(Ephys(34).sub_traces_train(:,2)))*1.1],'-','Color','b','LineWidth',2);
  hold on;
  end
%hold on;text(-550*srF,Ephys(28).IV.RMP,[num2str(rmp(28)),'mV'],'FontSize',9);
hold on;title('n-CPN');
subplot(1,3,3)
plot(Ephys(86).sub_traces_train(1:100000,2),'Color','m','LineWidth',1);set(gca,'box','off');axis off;
hold on;%ylim([ov_min-10 ov_max]);
  for h=1:5
  x1=delay+duration/freq*co(h);
  x2=(delay+pulsedur)+duration/freq*co(h);
  p1=plot([x1 x2],[abs(max(Ephys(86).sub_traces_train(:,2)))*1.1 abs(max(Ephys(86).sub_traces_train(:,2)))*1.1],'-','Color','b','LineWidth',2);
  hold on;
  end
% hold on;text(-550*srF,Ephys(15).IV.RMP,[num2str(rmp(15)),'mV'],'FontSize',9);
hold on;title('n-CPN IN');

%% Representative examples pA 
fig5=figure;set(fig5, 'Position', [200, 800, 600, 200]);set(gcf,'color','w');
subplot(1,3,1)
plot(Ephys(11).sub_traces_train(1:100000,1),'Color','g','LineWidth',1);set(gca,'box','off');axis off;
hold on;ylim([-230 30]);
freq=1;
co=[0:1:5-1];
duration=20000;
pulsedur=1000;
delay=5000;
  for h=1:5
  x1=delay+duration/freq*co(h);
  x2=(delay+pulsedur)+duration/freq*co(h);
  p1=plot([x1 x2],[abs(max(Ephys(11).sub_traces_train(:,2)))*3 abs(max(Ephys(11).sub_traces_train(:,2)))*3],'-','Color','b','LineWidth',2);
  hold on;
  end
%hold on;text(-550*srF,Ephys(7).IV.RMP,[num2str(rmp(7)),'mV'],'FontSize',9);
hold on;title('CPN');
subplot(1,3,2)
plot(Ephys(51).sub_traces_train(1:100000,1),'Color','k','LineWidth',1);set(gca,'box','off');axis off;
hold on;ylim([-230 30])
  for h=1:5
  x1=delay+duration/freq*co(h);
  x2=(delay+pulsedur)+duration/freq*co(h);
  p1=plot([x1 x2],[abs(max(Ephys(51).sub_traces_train(:,2)))*3 abs(max(Ephys(51).sub_traces_train(:,2)))*3],'-','Color','b','LineWidth',2);
  hold on;
  end
%hold on;text(-550*srF,Ephys(28).IV.RMP,[num2str(rmp(28)),'mV'],'FontSize',9);
hold on;title('n-CPN');
subplot(1,3,3)
plot(Ephys(37).sub_traces_train(1:100000,1),'Color','m','LineWidth',1);set(gca,'box','off');axis off;
hold on;%ylim([ov_min-10 ov_max]);
  for h=1:5
  x1=delay+duration/freq*co(h);
  x2=(delay+pulsedur)+duration/freq*co(h);
  p1=plot([x1 x2],[abs(max(Ephys(37).sub_traces_train(:,2)))*3 abs(max(Ephys(37).sub_traces_train(:,2)))*3],'-','Color','b','LineWidth',2);
  hold on;
  end

hold on;title('n-CPN IN');
ylim([-230 30]);

%% Quantification of mV and pA
%par=rmp;
% yl='Evoked Input (mV)';
% %yl='Rheobase (pA)';
% comi=[nonzeros(la_vc); nonzeros(nla_vc); nonzeros(in_vc)];
% group=[ones(1,length(nonzeros(la_vc))) ones(1,length(nonzeros(nla_vc)))*2 ones(1,length(nonzeros(in_vc)))*3];
% figure; set(gcf,'color','w');
% plot(group(1:length(nonzeros(la_vc))),nonzeros(la_vc),'go');hold on;plot(ones(1,length(nonzeros(nla_vc)))*2,nonzeros(nla_vc),'ko');hold on;plot(ones(1,length(nonzeros(in_vc)))*3,nonzeros(in_vc),'mo')
% %plotSpread(comi'-90,'categoryIdx',group, 'categoryMarkers',{'.','.'},'categoryColors',{'k','m'});
% hold on;
% boxplot(comi,group,'Colors','k','OutlierSize',0.001);
% box off;
% xticklabels({'CPNs','n-CPNs','n-CPN IN'});ylabel(yl);
figure; set(gcf,'color','w');
b1=bar(1,abs(nanmean(nonzeros(la_vc))));hold on;b2=bar(2.5,abs(nanmean(nonzeros(nla_vc))));b3=bar(4,abs(nanmean(nonzeros(in_vc))));
box off;xticks([1:1.5:4]);set(gca,'XTickLabel',{'CPNs' 'n-CPNs','n-CPN IN'});hold on;
b1.FaceColor=[0 1 0];b2.FaceColor=[0.3 0.3 0.3];b3.FaceColor=[0.8 0.2 1];
plot(ones(1,length(nonzeros(la_vc))),abs(nonzeros(la_vc)),'ko','MarkerEdgeColor',[0.7,0.7,0.7]);hold on;plot(ones(1,length(nonzeros(nla_vc)))*2.5,abs(nonzeros(nla_vc)),'ko','MarkerEdgeColor',[0.7,0.7,0.7]);hold on;plot(ones(1,length(nonzeros(in_vc)))*4,abs(nonzeros(in_vc)),'ko','MarkerEdgeColor',[0.7,0.7,0.7])
ylabel('Synaptic input (pA)');
err=nanstd(nonzeros(la_vc))/sqrt(length(nonzeros(la_vc)));
er = errorbar(1,abs(nanmean(nonzeros(la_vc))),err);    
er.Color = [0 0 0];      
er.LineWidth=1.5
er.LineStyle = 'none'; hold on;
err=nanstd(nonzeros(nla_vc))/sqrt(length(nonzeros(nla_vc)));
er = errorbar(2.5,abs(nanmean(nonzeros(nla_vc))),err);    
er.Color = [0 0 0]; 
er.LineWidth=1.5
er.LineStyle = 'none'; hold on;
err=nanstd(nonzeros(in_vc))/sqrt(length(nonzeros(in_vc)));
er = errorbar(4,abs(nanmean(nonzeros(in_vc))),err);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';hold on;
er.LineWidth=1.5;
%xtickangle(45);
ylim([0 580]);
yticks([50:100:575]);
hold on;
FS_spik=find(in_cc>30);
plot(ones(1,length(nonzeros(in_vc(FS_spik))))*4,abs(nonzeros(in_vc(FS_spik))),'ro','MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,0,0])
%% 
%% Quantification of time to peak
%par=rmp;
% yl='Evoked Input (mV)';
% %yl='Rheobase (pA)';
% comi=[nonzeros(la_vc); nonzeros(nla_vc); nonzeros(in_vc)];
% group=[ones(1,length(nonzeros(la_vc))) ones(1,length(nonzeros(nla_vc)))*2 ones(1,length(nonzeros(in_vc)))*3];
% figure; set(gcf,'color','w');
% plot(group(1:length(nonzeros(la_vc))),nonzeros(la_vc),'go');hold on;plot(ones(1,length(nonzeros(nla_vc)))*2,nonzeros(nla_vc),'ko');hold on;plot(ones(1,length(nonzeros(in_vc)))*3,nonzeros(in_vc),'mo')
% %plotSpread(comi'-90,'categoryIdx',group, 'categoryMarkers',{'.','.'},'categoryColors',{'k','m'});
% hold on;
% boxplot(comi,group,'Colors','k','OutlierSize',0.001);
% box off;
% xticklabels({'CPNs','n-CPNs','n-CPN IN'});ylabel(yl);
srf=20;
figure; set(gcf,'color','w');
b1=bar(1,abs(nanmean(nonzeros(la_vc_d)/srf)));hold on;b2=bar(2.5,abs(nanmean(nonzeros(nla_vc_d)/srf)));b3=bar(4,abs(nanmean(nonzeros(in_vc_d)/srf)));
box off;xticks([1:1.5:4]);set(gca,'XTickLabel',{'CPNs' 'n-CPNs','n-CPN IN'});hold on;
b1.FaceColor=[0 1 0];b2.FaceColor=[0.3 0.3 0.3];b3.FaceColor=[0.8 0.2 1];
plot(ones(1,length(nonzeros(la_vc_d))),abs(nonzeros(la_vc_d)/srf),'ko','MarkerEdgeColor',[0.7,0.7,0.7]);hold on;plot(ones(1,length(nonzeros(nla_vc_d)))*2.5,abs(nonzeros(nla_vc_d)/srf),'ko','MarkerEdgeColor',[0.7,0.7,0.7]);hold on;plot(ones(1,length(nonzeros(in_vc_d)))*4,abs(nonzeros(in_vc_d)/srf),'ko','MarkerEdgeColor',[0.7,0.7,0.7])
ylabel('Time to peak (ms)');
err=nanstd(nonzeros(la_vc_d)/srf)/sqrt(length(nonzeros(la_vc_d)/srf));
er = errorbar(1,abs(nanmean(nonzeros(la_vc_d)/srf)),err);    
er.Color = [0 0 0];      
er.LineWidth=1.5
er.LineStyle = 'none'; hold on;
err=nanstd(nonzeros(nla_vc_d)/srf)/sqrt(length(nonzeros(nla_vc_d)));
er = errorbar(2.5,abs(nanmean(nonzeros(nla_vc_d)/srf)),err);    
er.Color = [0 0 0]; 
er.LineWidth=1.5
er.LineStyle = 'none'; hold on;
err=nanstd(nonzeros(in_vc_d)/srf)/sqrt(length(nonzeros(in_vc_d)/srf));
er = errorbar(4,abs(nanmean(nonzeros(in_vc_d)/srf)),err);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';hold on;
er.LineWidth=1.5;
%xtickangle(45);
ylim([0 50]);
yticks([0:10:50]);
%% 

 for i=1:length(nla)
     if Ephys(nla(i)).wc==1 & Ephys(nla(i)).sol==1
     if size(Ephys(nla(i)).train_p,2)>=2
    nla_cc_n(i)=min(Ephys(nla(i)).train_n(:,2));
     else
     nla_cc_n(i)=NaN;
     end
     else
     nla_cc_n(i)=NaN;
     end
 end
 %% Look at Cs-gluc
 for i=1:size(Ephys,2)
 sol(i)=Ephys(i).sol;
 end
cs=find(sol==2);
%% 
fig5=figure;set(fig5, 'Position', [200, 800, 400, 200]);set(gcf,'color','w');
plot(Ephys(cs(3)).sub_traces_train(1:100000,1),'Color','r','LineWidth',1);set(gca,'box','off');axis off;
hold on;plot(Ephys(cs(3)).sub_traces_train(1:100000,2)+20,'Color','b','LineWidth',1);set(gca,'box','off');axis off;
hold on;
freq=1;
co=[0:1:5-1];
duration=20000;
pulsedur=1000;
delay=5000;
  for h=1:5
  x1=delay+duration/freq*co(h);
  x2=(delay+pulsedur)+duration/freq*co(h);
  p1=plot([x1 x2],[abs(max(Ephys(11).sub_traces_train(:,2)))*10 abs(max(Ephys(11).sub_traces_train(:,2)))*10],'-','Color','b','LineWidth',2);
  hold on;
  end
%hold on;text(-550*srF,Ephys(7).IV.RMP,[num2str(rmp(7)),'mV'],'FontSize',9);

%% 
 %% Look at Cs-gluc
 for i=1:size(Ephys,2)
 tmp(i)=Ephys(i).wc;
 end
wc=find(tmp==0);

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