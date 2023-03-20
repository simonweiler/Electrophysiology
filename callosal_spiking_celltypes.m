%% Load data structure 
str_L6    = 'D:\Postdoc_Margrie\Projects\Callosal\output';
folder_list = uipickfiles('FilterSpec',str_L6);
load(char(folder_list));
%sampling rate
srF=20;
sr=20000;
%% Load GAD, PV and SOM with K-gluconate
%GAD
gad=cell_selecter(Ephys,'drugs',0,'label',2,'geno',4,'sol',1);
%PV
pv=cell_selecter(Ephys,'drugs',0,'label',3,'geno',5,'sol',1);
%SOM
som=cell_selecter(Ephys,'drugs',0,'label',6,'geno',8,'sol',1);
%% Read out light evoked epsp amplitude (including spikes)
[epsp_gad] = readout_amp_epsp(Ephys,gad,1,sr);
[epsp_pv] = readout_amp_epsp(Ephys,pv,2,sr);
[epsp_som] = readout_amp_epsp(Ephys,som,2,sr);
%% Read out max spike frequency
 temp=[];maxsF_gad=[];
temp=find(gad==1);
spike_counts_all=[];
for i=1:length(temp)
    try
    maxsF_gad(i)=max(Ephys(temp(i)).IV.spikecount);
    spikecount_gad{:,i}=Ephys(temp(i)).IV.spikecount;
    temp3=Ephys(temp(i)).IV.spikecount;
    spike_countsm=ones(1,25)*NaN;
    if length(temp3)==15
    spike_countsm(1:15)=temp3;
    else
    spike_countsm(1:25)=temp3;
    end
    stimvec_gad{:,i}=Ephys(temp(i)).IV.stimvec;
    spike_counts_all(:,i)=spike_countsm;

    catch
    maxsF_gad(i)=NaN;
    spikecount_gad{:,i}=NaN;
    stimvec_gad{:,i}=NaN;
    spike_counts_all(:,i)=ones(1,25)*NaN;
    end
end
%% 
 temp=[];maxsF_pv=[];
temp=find(pv==1);
spike_counts_all=[];
for i=1:length(temp)
    try
    maxsF_pv(i)=max(Ephys(temp(i)).IV.spikecount);
    spikecount_pv{:,i}=Ephys(temp(i)).IV.spikecount;
    temp3=Ephys(temp(i)).IV.spikecount;
    spike_countsm=ones(1,25)*NaN;
    if length(temp3)==15
    spike_countsm(1:15)=temp3;
    else
    spike_countsm(1:25)=temp3;
    end
    stimvec_pv{:,i}=Ephys(temp(i)).IV.stimvec;
    spike_counts_all(:,i)=spike_countsm;

    catch
    maxsF_pv(i)=NaN;
    spikecount_pv{:,i}=NaN;
    stimvec_pv{:,i}=NaN;
    spike_counts_all(:,i)=ones(1,25)*NaN;
    end
end
%% 
 temp=[];maxsF_som=[];
temp=find(som==1);
spike_counts_all=[];
for i=1:length(temp)
    try
    maxsF_som(i)=max(Ephys(temp(i)).IV.spikecount);
    spikecount_som{:,i}=Ephys(temp(i)).IV.spikecount;
    temp3=Ephys(temp(i)).IV.spikecount;
    spike_countsm=ones(1,25)*NaN;
    if length(temp3)==15
    spike_countsm(1:15)=temp3;
    else
    spike_countsm(1:25)=temp3;
    end
    stimvec_som{:,i}=Ephys(temp(i)).IV.stimvec;
    spike_counts_all(:,i)=spike_countsm;

    catch
    maxsF_som(i)=NaN;
    spikecount_som{:,i}=NaN;
    stimvec_som{:,i}=NaN;
    spike_counts_all(:,i)=ones(1,25)*NaN;
    end
end
%% 
[rmp_gad maxsp_gad rheo_gad rin_gad tau_gad sag_gad trace_gad spike_time_gad] = passive_readout(Ephys,gad);
[rmp_pv maxsp_pv rheo_pv rin_pv tau_pv sag_pv trace_pv spike_time_pv] = passive_readout(Ephys,pv);
[rmp_som maxsp_som rheo_som rin_som tau_som sag_som trace_som spike_time_som] = passive_readout(Ephys,som);


%% 
spike_nr=1;
active_gad=sp_parameters_pandora(trace_gad,spike_nr);
active_pv=sp_parameters_pandora(trace_pv,spike_nr);
active_som=sp_parameters_pandora(trace_som,spike_nr);
%% 
nrt=6;
comp_gad=active_gad(nrt,:);
comp_pv=active_pv(nrt,:);
comp_som=active_som(nrt,:);
figure;scatter(comp_som,epsp_som);hold on;scatter(comp_pv,epsp_pv);hold on;scatter(comp_gad,epsp_gad)