%Callosal input mapping and ephys charaterization  scripts 
%genotypes 1=Tlx3Cre; 2=Rpb4Cre; 3=Syt6Cre; 4=GADcre tdtomato; 5= PVcre tdtomat; 6=Ntsr1cre tdtomato 
%% Load data structure 
str_L6    = 'C:\Users\slice setup\Electrophysiology\output';
folder_list = uipickfiles('FilterSpec',str_L6);
load(char(folder_list));
%sampling rate
srF=20;
sr=20000;
%% use filter function 'cell selecter' to read out desired cells/line etc.
%Ntsr1 mouse line, K-gluc, retro cells
retro_k = cell_selecter(Ephys, 'label',1, 'sol',1,'geno',6);
%Ntsr1 mouse line, K-gluc, ntsr1+ cells
ntsr_k = cell_selecter(Ephys, 'label',4, 'sol',1,'geno',6);
%Ntsr1 mouse line, Cs-gluc, retro cells
retro_cs = cell_selecter(Ephys, 'label',1, 'sol',2,'geno',6);
%Ntsr1 mouse line, Cs-gluc, ntsr1+ cells
ntsr_cs = cell_selecter(Ephys, 'label',4, 'sol',2,'geno',6);

%PV mouse line, K-gluc, retro cells
retro_k_pv = cell_selecter(Ephys, 'label',1, 'sol',1,'geno',5);
%PV mouse line, K-gluc, ntsr1+ cells
pv_k = cell_selecter(Ephys, 'label',3, 'sol',1,'geno',5);
%PV mouse line, Cs-gluc, retro cells
retro_cs_pv = cell_selecter(Ephys, 'label',1, 'sol',2,'geno',5);
%PV mouse line, Cs-gluc, ntsr1+ cells
pv_cs = cell_selecter(Ephys, 'label',3, 'sol',2,'geno',5);
%% use filter for all retro cells with K-gluc
rk = cell_selecter(Ephys, 'label',1, 'sol',1);
%% use filter for all Interneuorns with K-gluc
gadk = cell_selecter(Ephys, 'label',[2], 'sol',1);
%% 
nlk = cell_selecter(Ephys, 'label',[0], 'sol',1);
%% 100% retro are ntrs-
fig3= figure;set(fig3, 'Name', 'Paired comp');set(fig3, 'Position', [200, 300, 250, 200]);set(gcf,'color','w');
fdat=[100 0];
b=bar([1 2],fdat);ylabel('Fraction CPN (%)');xticklabels({'Ntsr1(-)', 'Ntsr1(+)'});box off;
b.FaceColor=[0.5 0.5 0.5];set(gca,'FontSize',10);
%% Plot IV for all retro cells
plot_intrinsic(Ephys,rk,8,srF,'k',1);
%% Plot IV for all GAD+ cells
plot_intrinsic(Ephys,gadk,6,srF,'m',1);
%% Plot IV for all PV+ cells
plot_intrinsic(Ephys,pv_k,2,srF,'m',1);
%% Plot IV for all Ntsr1+ cells
plot_intrinsic(Ephys,ntsr_k,4,srF,'r',1);
%% Plot Rheobase 1x and 2x for all Ntsr1+ cells
[F2xRheo]=plot_intrinsic(Ephys,ntsr_k,4,srF,'k',2);
%% Plot Rheobase 1x and 2x for all retro+ cells
[F2xRheo_rk]=plot_intrinsic(Ephys,retro_k,4,srF,'m',2);
%% Plot Rheobase 1x and 2x for all nonlabelled cells
plot_intrinsic(Ephys,nlk,5,srF,'m',2);
%% 
tim_1=[F2xRheo_rk' F2xRheo'];

%% Instrinisc properties Ntsr1 vs retro
%extract info from IV trace and passive, extract Rheobase traces where at
%least 2 spikes are present
%retro
temp=[];
for i=1:length(find(retro_k==1));
    temp=find(retro_k==1);
    rmp_retro(i)=Ephys(temp(i)).IV.RMP;
    maxsp_retro(i)=max(Ephys(temp(i)).IV.spikecount);
    rheo_retro(i)=Ephys(temp(i)).Rheobase.rheo;
    rin_retro(i)=Ephys(temp(i)).Passive.Rin;
    tau_retro(i)=Ephys(temp(i)).Passive.tau(1);
    if isempty(Ephys(temp(i)).Sag)==0
    sag_retro(i)=Ephys(temp(i)).Sag.sagratio(1);
    else
    sag_retro(i)=NaN;
    end
    md=[];
     if isempty(find(Ephys(temp(i)).Rheobase.spikecount>=2))==0
    md=find(Ephys(temp(i)).Rheobase.spikecount>=2);
    trace_retro(:,i)=Ephys(temp(i)).Rheobase.traces(:,md(1));
   else
      md=find(Ephys(temp(i)).Rheobase.spikecount>=1);
    trace_retro(:,i)=Ephys(temp(i)).Rheobase.traces(:,md(1));
   end
end
%Ntsr1
temp=[];
for i=1:length(find(ntsr_k==1));
    temp=find(ntsr_k==1);
    rmp_ntsr(i)=Ephys(temp(i)).IV.RMP;
    maxsp_ntsr(i)=max(Ephys(temp(i)).IV.spikecount);
    rheo_ntsr(i)=Ephys(temp(i)).Rheobase.rheo;
    rin_ntsr(i)=Ephys(temp(i)).Passive.Rin;
    tau_ntsr(i)=Ephys(temp(i)).Passive.tau(1);
    if isempty(Ephys(temp(i)).Sag)==0
    sag_ntsr(i)=Ephys(temp(i)).Sag.sagratio(1);
    else
    sag_ntsr(i)=NaN;
    end
    md=[];
   if isempty(find(Ephys(temp(i)).Rheobase.spikecount>=2))==0
    md=find(Ephys(temp(i)).Rheobase.spikecount>=2);
    trace_ntsr(:,i)=Ephys(temp(i)).Rheobase.traces(:,md(1));
   else
      md=find(Ephys(temp(i)).Rheobase.spikecount>=1);
    trace_ntsr(:,i)=Ephys(temp(i)).Rheobase.traces(:,md(1));
   end
end
rmp_1= [rmp_retro' rmp_ntsr'];
maxsp_1=[maxsp_retro' maxsp_ntsr'];
rheo_1=[rheo_retro' rheo_ntsr'];
rin_1=[rin_retro' rin_ntsr'];
tau_1=[tau_retro' tau_ntsr'];
sag_1=[sag_retro' sag_ntsr'];
%% Reading out further active spike properties using Pandora
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
%Pandora parameters
sp=2;
temp=[];
for i=1:length(find(retro_k==1));
trace_curr=trace_retro(:,i);
dt = 1e-4;dy = 1e-3;
props = struct('spike_finder', 2, 'threshold', 0);
traces_analysis=trace(trace_curr,dt, dy, 'Analysis', props);
alltrace_info = getProfileAllSpikes(traces_analysis);
parameters=[];
parameters=alltrace_info.spikes_db.data; 
try
active_ntsr(:,i)=parameters(sp,:);
catch
active_ntsr(:,i)=parameters(1,:);
end;
end
temp=[];
for i=1:length(find(ntsr_k==1));
trace_curr=trace_ntsr(:,i);
dt = 1e-4;dy = 1e-3;
props = struct('spike_finder', 2, 'threshold', 0);
traces_analysis=trace(trace_curr,dt, dy, 'Analysis', props);
alltrace_info = getProfileAllSpikes(traces_analysis);
parameters=[];
parameters=alltrace_info.spikes_db.data;
try
active_ntsr(:,i)=parameters(sp,:);
catch
active_ntsr(:,i)=parameters(1,:);
end
end
%% Test rmp, etc
close all;
%RMP
paired_plot(rmp_1,0);xticklabels({'CPN','NtsR1'});ylabel('RMP (mV)');yticks([-80:5:-60]);set(gca,'FontSize',10);
%Input resistance Rin
paired_plot(rin_1,0);xticklabels({'CPN','NtsR1'});ylabel('Input resistance (mhm)');set(gca,'FontSize',10);
%Tau 
paired_plot(tau_1,0);xticklabels({'CPN','NtsR1'});ylabel('Tau (ms)');set(gca,'FontSize',10);
%Maximum spike rate
paired_plot(maxsp_1,0);xticklabels({'CPN','NtsR1'});ylabel('Max spike number');yticks([0:5:30]);set(gca,'FontSize',10);
%Rheobase
paired_plot(rheo_1,1);xticklabels({'CPN','NtsR1'});ylabel('Rheobase (pA)');yticks([0:50:200]);set(gca,'FontSize',10);
%Sag Ratio
paired_plot(sag_1,1);xticklabels({'CPN','NtsR1'});ylabel('Sag Ratio');set(gca,'FontSize',10);
%% Pired timing of first two spikes

paired_plot(tim_1,1);xticklabels({'CPN','NtsR1'});ylabel('initial spikes');set(gca,'FontSize',10);
%% Plotting paired using the second spike always (if possible) 
close all;
paired_subplot(active_retro([1 2 3 4 5 6 7 8 10 11 12 14 15],:),active_ntsr([1 2 3 4 5 6 7 8 10 11 12 14 15],:),1)
%% Input CPN vs Ntsr1
temp=[];
for i=1:length(find(retro_cs==1));
    temp=find(retro_cs==1);
    retro_epsc(i)=max(abs(Ephys(temp(i)).train_n(1,:)));
    retro_ipsc(i)=max(abs(Ephys(temp(i)).train_p(1,:)));
end
temp=[];
for i=1:length(find(ntsr_cs==1));
    temp=find(ntsr_cs==1);
    ntsr_epsc(i)=max(abs(Ephys(temp(i)).train_n(1,:)));
    ntsr_ipsc(i)=max(abs(Ephys(temp(i)).train_p(1,:)));
end
temp=[];
for i=1:length(find(retro_k==1));
    temp=find(retro_k==1);
    retro_epsc2(i)=max(abs(Ephys(temp(i)).train_n(1,:)));    
end
temp=[];
for i=1:length(find(ntsr_k==1));
    temp=find(ntsr_k==1);
    ntsr_epsc2(i)=max(abs(Ephys(temp(i)).train_n(1,:)));    
end
epsc_1=[[retro_epsc retro_epsc2(9:12)]' [ntsr_epsc ntsr_epsc2(9:12)]'];
ipsc_1=[retro_ipsc' ntsr_ipsc'];
e_i_retro=retro_epsc./retro_ipsc;
e_i_ntsr=ntsr_epsc./ntsr_ipsc;
e_i_1=[e_i_retro' e_i_ntsr'];
%% Plot EPSC and IPSC retro vs NtsR1
%EPSC
paired_plot(epsc_1,0);xticklabels({'CPN','NtsR1'});ylabel('EPSC peak (pA)');set(gca,'FontSize',10);
%IPSC
paired_plot(ipsc_1,0);xticklabels({'CPN','NtsR1'});ylabel('IPSC peak (pA)');set(gca,'FontSize',10);
%E_I
paired_plot(e_i_1,0);xticklabels({'CPN','NtsR1'});ylabel('E / I Ratio');set(gca,'FontSize',10);
%% 



%% PV cell vs retro
temp=[];
for i=1:length(find(retro_cs_pv==1));
    temp=find(retro_cs_pv==1);
    retro_epsc_pv(i)=max(abs(Ephys(temp(i)).train_n(1,:)));
    retro_ipsc_pv(i)=max(abs(Ephys(temp(i)).train_p(1,:)));
end
temp=[];
for i=1:length(find(pv_cs==1));
    temp=find(pv_cs==1);
    pv_epsc(i)=max(abs(Ephys(temp(i)).train_n(1,:)));
    pv_ipsc(i)=max(abs(Ephys(temp(i)).train_p(1,:)));
end
temp=[];
for i=1:length(find(retro_k_pv==1));
    temp=find(retro_k_pv==1);
    retro_epsc2_pv(i)=max(abs(Ephys(temp(i)).train_n(1,:)));    
end
temp=[];
for i=1:length(find(pv_k==1));
    temp=find(pv_k==1);
    pv_epsc2(i)=max(abs(Ephys(temp(i)).train_n(1,:)));    
end
epsc_2=[[retro_epsc_pv retro_epsc2_pv]' [pv_epsc pv_epsc2]'];
ipsc_2=[retro_ipsc_pv' pv_ipsc'];
e_i_retro_pv=retro_epsc_pv./retro_ipsc_pv;
e_i_pv=pv_epsc./pv_ipsc;
e_i_2=[e_i_retro_pv' e_i_pv'];
%% Plot EPSC and IPSC retro vs PV
%EPSC
paired_plot(epsc_2,0);xticklabels({'CPN','PV'});ylabel('EPSC peak (pA)');set(gca,'FontSize',10);
%IPSC
paired_plot(ipsc_2,0);xticklabels({'CPN','PV'});ylabel('IPSC peak (pA)');set(gca,'FontSize',10);
%E_I
paired_plot(e_i_2,0);xticklabels({'CPN','PV'});ylabel('E / I Ratio');set(gca,'FontSize',10);
