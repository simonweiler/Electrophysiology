%Callosal input mapping and ephys charaterization  scripts 
%genotypes 1=Tlx3Cre; 2=Rpb4Cre; 3=Syt6Cre; 4=GADcre tdtomato; 5= PVcre
%tdtomat; 6=Ntsr1cre tdtomato ;7=Penk; 8=SOM
%labels
%1=callosal; 2=GAD; 3=PV; 4=Ntsr1 tdtomato; 5= Penk
%6=SOM
%% Load data structure 
str_L6    = 'D:\Postdoc_Margrie\Projects\Callosal\output';
folder_list = uipickfiles('FilterSpec',str_L6);
load(char(folder_list));
%sampling rate
srF=20;
sr=20000;
%% read out SOM cells with Cs in SOM animals 
som_cs = cell_selecter(Ephys, 'label',6,'geno',8,'sol',2,'optovariant',1);
%% read out CPNs cells with Cs in SOM animals
cpn_som_cs = cell_selecter(Ephys, 'label',1,'geno',8,'sol',2,'optovariant',1);
%% 
som_k = cell_selecter(Ephys, 'label',6,'geno',8,'sol',1,'optovariant',1);
%% 
cpn_som_k = cell_selecter(Ephys, 'label',1,'geno',8,'sol',1,'optovariant',1);
%% read out EPSCs/IPSCSs SOM Cs in SOM animals
[esom_cs_long isom_cs_long eisom_cs_long] = readout_amp(Ephys,som_cs ,1,2,1,2);
%% read out EPSCs/IPSCSs CPNs Cs in SOM animals 
[ecpnsom_cs_long icpnsom_cs_long eicpnsom_cs_long] = readout_amp(Ephys,cpn_som_cs ,1,2,1,2);
%% read out EPSCs/IPSCSs SOM Cs in SOM animals
[esom_k_long isom_k_long eisom_k_long] = readout_amp(Ephys,som_k ,1,1,1,2);
%% read out EPSCs/IPSCSs CPNs Cs in SOM animals 
[ecpnsom_k_long icpnsom_k_long eicpnsom_k_long] = readout_amp(Ephys,cpn_som_k ,1,1,1,2);
%% 

temp1=[];temp2=[];
for i=1:6
temp1(i,:) = cell_selecter(Ephys,'label',[4],'sol',2,'geno',6,'optovariant',1,'pair',i,'drugs',1);
temp2(i,:) = cell_selecter(Ephys,'label',[1],'sol',2,'geno',6,'optovariant',1,'pair',i,'drugs',1);
end
cre_on_cs=sum(temp1);cre_off_cs=sum(temp2);
%% read out EPSCs IPSCs in Ntsre Cs
ntsr_cs = cell_selecter(Ephys, 'label',4,'geno',6,'sol',2,'optovariant',1,'drugs',0);
ntsr_cs_wi = cell_selecter(Ephys, 'label',4,'geno',6,'sol',2,'optovariant',1,'drugs',1);
%% 
cpns_ntsr_cs = cell_selecter(Ephys, 'label',1,'geno',6,'sol',2,'optovariant',1,'drugs',0);
%cpns_ntsr_cs_wi = cell_selecter(Ephys, 'label',1,'geno',6,'sol',2,'optovariant',1,'drugs',1);
%% 

[entsr_cs_hf intsr_cs_hf eintsr_cs_hf] = readout_amp(Ephys,ntsr_cs ,1,2,1,2);
[entsr_cs_hf_w intsr_cs_hf_w eintsr_cs_hf_w] = readout_amp(Ephys,ntsr_cs_wi ,1,2,1,2);
%% 
[ecpnr_ntsr_cs_hf icpnr_ntsr_cs_hf eicpnr_ntsr_cs_hf] = readout_amp(Ephys,cpns_ntsr_cs ,1,2,1,2);
%[ecpnr_ntsr_cs_hf_w icpnr_ntsr_cs_hf_w eicpnr_ntsr_cs_hf_w] = readout_amp(Ephys,cpns_ntsr_cs_wi ,1,2,1,2);
%% 
com_cs=[];
com_cs=[[eintsr_cs_hf eintsr_cs_hf_w(1:2)]' eicpnr_ntsr_cs_hf'];
[gh gm]=find(isnan(com_cs));
com_cs(unique(gh),:)=[];
cl={'r','b'};
data=[];data=com_cs;
paired_plot_box(data,cl);
xticklabels({'EX','IN'});ylabel('Onset Latency (ms)');set(gca,'FontSize',10);
%% 

com_cpnsom_cs=[isom_cs_long' icpnsom_cs_long'];
[gh gm]=find(isnan(com_cpnsom_cs));
com_cpnsom_cs(unique(gh),:)=[];
cl={'r','b'};
data=[];data=com_cpnsom_cs;
paired_plot_box(data,cl);
xticklabels({'EX','IN'});ylabel('Onset Latency (ms)');set(gca,'FontSize',10);
%% 
com_som=[]
com_som=[[esom_cs_long esom_k_long]' [ecpnsom_cs_long ecpnsom_k_long]'];
cl={'r','b'};
data=[];data=com_som;
paired_plot_box(data,cl);
xticklabels({'EX','IN'});ylabel('Onset Latency (ms)');set(gca,'FontSize',10);

%% use filter function 'cell selecter' to read out desired cells/line etc.
%Ntsr1 mouse line, K-gluc, retro cells
%retro_k = cell_selecter(Ephys, 'label',1, 'sol',1,'geno',6);
%Ntsr1 mouse line, K-gluc, ntsr1+ cells
ntsr_opto1_cs = cell_selecter(Ephys, 'label',4,'geno',6,'sol',2,'optovariant',0);
%% 
ntsr_all_cs = cell_selecter(Ephys, 'label',4,'geno',6,'sol',2);
%% 
retro_all_cs = cell_selecter(Ephys, 'label',1,'geno',6,'sol',2);

%% 

ntsr_ttx_k = cell_selecter(Ephys, 'label',4,'geno',6, 'drugs',2,'sol',1);
ntsr_ttx_cs = cell_selecter(Ephys, 'label',4,'geno',6, 'drugs',2,'sol',2);
ntsr_wttx_k= cell_selecter(Ephys, 'label',4,'geno',6, 'drugs',1,'sol',1);
ntsr_wttx_cs= cell_selecter(Ephys, 'label',4,'geno',6, 'drugs',1,'sol',2);
%% penk mouse line 
penk_cs = cell_selecter(Ephys, 'label',5,'geno',7,'sol',2,'optovariant',1);
retro_cs_penk = cell_selecter(Ephys, 'label',1,'geno',7,'sol',2,'optovariant',1);
%% 
all_retro_cs = cell_selecter(Ephys, 'label',1,'sol',2,'drugs',0);
%% 
all_ntsr1_cs = cell_selecter(Ephys, 'label',4,'sol',2,'optovariant',1,'drugs',1);



%% 
[epsc_all_retro_long ipsc_all_retro_long e_i_all_retro_train] = readout_amp(Ephys,all_retro_cs ,1,2,1,2);
%% 
[epsc_all_ntsr1_long ipsc_all_ntsr1_long e_i_all_ntsr1_train] = readout_amp(Ephys,all_ntsr1_cs ,1,2,1,2);

%% 
[epsc_retropenk_long ipsc_retropenk_long e_i_retropenk_train] = readout_amp(Ephys,retro_cs_penk ,1,2,1,2);
[epsc_retropenk_hf ipsc_retropenk_hf e_i_ratio_retropenk_hf] = readout_amp(Ephys,retro_cs_penk ,2,2,1,2);
[epsc_retropenk_hf2 ipsc_retropenk_hf2 e_i_ratio_retropenk_hf2] = readout_amp(Ephys,retro_cs_penk,3,2,1,2);
%% 
[epsc_penk_long ipsc_penk_long e_i_penk_train] = readout_amp(Ephys,penk_cs ,1,2,1,2);
[epsc_penk_hf ipsc_penk_hf e_i_ratio_penk_hf] = readout_amp(Ephys,penk_cs ,2,2,1,2);
[epsc_penk_hf2 ipsc_penk_hf2 e_i_ratio_penk_hf2] = readout_amp(Ephys,penk_cs,3,2,1,2);
%% 
ei_penk_retro=[epsc_penk_long' epsc_retropenk_long'];
[gh gm]=find(isnan(ei_penk_retro));
ei_penk_retro(unique(gh),:)=[];
cl={'r','b'};
data=[];data=ei_penk_retro;
paired_plot_box(data,cl);
xticklabels({'EX','IN'});ylabel('Onset Latency (ms)');set(gca,'FontSize',10);


%% 


%Ntsr1 mouse line, Cs-gluc, retro cells
%retro_cs = cell_selecter(Ephys, 'label',1, 'sol',2,'geno',6);
%Ntsr1 mouse line, Cs-gluc, ntsr1+ cells
%ntsr_cs = cell_selecter(Ephys, 'label',4, 'sol',2,'geno',6);
%% 
temp=[];cells_idx=[];
cnr=1;
cells_idx=ntsr_all_cs;
temp=find(cells_idx==1);
fig4=figure;set(fig4, 'Position', [200, 200, 300, 300]);set(gcf,'color','w');
plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,1),'Color','k','LineWidth',1);set(gca,'box','off');
hold on;plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,2),'Color','r','LineWidth',1);set(gca,'box','off');
%% 
temp=[];cells_idx=[];
cnr=1;
cells_idx=retro_all_cs;
temp=find(cells_idx==1);
fig4=figure;set(fig4, 'Position', [200, 200, 300, 300]);set(gcf,'color','w');
for i=1:length(temp)
plot(Ephys(temp(i)).sub_traces_train(1:1*sr,2),'LineWidth',1);set(gca,'box','off');
hold on;
end
%% 
temp=[];cells_idx=[];
cnr=4;
cells_idx=ntsr_wttx_cs;
temp=find(cells_idx==1);
fig4=figure;set(fig4, 'Position', [200, 200, 300, 300]);set(gcf,'color','w');
plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,1),'Color','k','LineWidth',1);set(gca,'box','off');
hold on;plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,3),'Color','r','LineWidth',1);set(gca,'box','off');
%% Read out max epsc and ipsc amplitude for train stimulus, and the two high frequency stimuli
cells_idx=[];cells_idx=ntsr_ttx;stim_type=1;
[train_e train_i train_ei] = readout_amp(Ephys,cells_idx,stim_type);
stim_type=2;
[hf_e hf_i hf_ei] = readout_amp(Ephys,cells_idx,stim_type);
stim_type=3;
[hf2_e hf2_i hf2_ei] = readout_amp(Ephys,cells_idx,stim_type);

% 
% %PV mouse line, K-gluc, retro cells
% retro_k_pv = cell_selecter(Ephys, 'label',1, 'sol',1,'geno',5);
% %PV mouse line, K-gluc, pv+ cells
% pv_k = cell_selecter(Ephys, 'label',3, 'sol',1,'geno',5);
%% 

%PV mouse line, Cs-gluc, retro cells
retro_cs_pv = cell_selecter(Ephys, 'label',1, 'sol',2,'geno',5);
%PV mouse line, Cs-gluc,PV+ cells
pv_cs = cell_selecter(Ephys, 'label',3, 'sol',2,'geno',5);
tb=[];ta=[];
tb=find(retro_cs_pv==1);
ta=find([Ephys(retro_cs_pv).pair]==0);
tb(ta);
retro_cs_pv(tb(ta))=0;

% tb=[];ta=[];
% tb=find(pv_k==1);
% ta=find([Ephys(pv_k).pair]==0);
% pv_k(tb(ta))=0;


%% 
% NtsR1 mouse line, K-gluc, non labelled
non_ntsr=cell_selecter(Ephys, 'label',0, 'sol',1,'geno',6);
% NtsR1 mouse line, K-gluc,  labelled but only paired with non ntsr
% unlabbeled
non_ntsr_lab=cell_selecter(Ephys, 'label',1, 'sol',1,'geno',6);
 te=find(non_ntsr_lab==1);
non_ntsr_lab(te(1:end-14))=0;
%% 
%% use filter for all retro cells with K-gluc
rk = cell_selecter(Ephys, 'label',1, 'sol',1);
%% use filter for all Interneuorns with K-gluc
ink1 = cell_selecter(Ephys, 'label',[2], 'sol',1);
ink2 = cell_selecter(Ephys, 'label',[3], 'sol',1);
ink=ink1+ink2;
%% nonlabelled cells
nlk = cell_selecter(Ephys, 'label',[0], 'sol',1);
%% retro cells in GAD mouse cs
retro_gadcs = cell_selecter(Ephys, 'label',[1], 'sol',2,'geno',4);
%% 
gadk = cell_selecter(Ephys, 'label',[2], 'sol',1,'geno',4);
%% 100% retro are ntrs-
fig3= figure;set(fig3, 'Name', 'Paired comp');set(fig3, 'Position', [200, 300, 250, 200]);set(gcf,'color','w');
fdat=[100 0];
b=bar([1 2],fdat);ylabel('Fraction CPN (%)');xticklabels({'Ntsr1(-)', 'Ntsr1(+)'});box off;
b.FaceColor=[0.5 0.5 0.5];set(gca,'FontSize',12);xtickangle(45)
%% 100% retro are GAD-
fig3= figure;set(fig3, 'Name', 'Paired comp');set(fig3, 'Position', [200, 300, 250, 200]);set(gcf,'color','w');
fdat=[100 0];
b=bar([1 2],fdat);ylabel('Fraction CPN (%)');xticklabels({'GAD(-)', 'GAD(+)'});box off;
b.FaceColor=[0.5 0.5 0.5];set(gca,'FontSize',12);xtickangle(45)
%% Plot IV for all retro cells
plot_intrinsic(Ephys,rk,8,srF,'k',1);
%% Plot IV for all GAD+ cells
plot_intrinsic(Ephys,gadk,6,srF,'m',1);
%% Plot IV for all PV+ cells
plot_intrinsic(Ephys,pv_k,2,srF,'m',1);
%% Plot IV for all Ntsr1+ cells
plot_intrinsic(Ephys,ntsr_k,4,srF,'r',1);
%% Plot IV for all reto in ntsr mouse line cells
plot_intrinsic(Ephys,retro_k,4,srF,'k',1);
%% Plot Rheobase 1x and 2x for all Ntsr1+ cells
[F2xRheo]=plot_intrinsic(Ephys,ntsr_k,4,srF,'k',2);
%% Plot Rheobase 1x and 2x for all retro+ cells
[F2xRheo_rk]=plot_intrinsic(Ephys,retro_k,4,srF,'m',2);
%% Plot Rheobase 1x and 2x for all nonlabelled cells
plot_intrinsic(Ephys,nlk,5,srF,'m',2);
%% Time fraction between first two spikes
tim_1=[F2xRheo_rk' F2xRheo'];
%% Read out intrinsic properties for all retro cells with K-gluc across all mouse lines
%extract info from IV trace and passive, extract Rheobase traces where at
%least 2 spikes are present
[rmp_rk maxsp_rk rheo_rk rin_rk tau_rk sag_rk trace_rk st_rk] = passive_readout(Ephys,rk);
active_rk=sp_parameters_pandora(trace_rk,2);
%% Plot histogram
rk_all=vertcat(rmp_rk,rin_rk,tau_rk,maxsp_rk,rheo_rk,sag_rk,active_rk([1 2 3 4 5 6 7 8 10 11 14 15],:));
fig9=figure;set(fig9, 'Position', [200, 200, 1400, 400]);set(gcf,'color','w');
str={'RMP (mV)','R_{IN}(M\Omega)','Tau (ms)','Max Nr. Spikes','Rheobase (pA)','Sag ratio','V_{min} (mV)','V_{peak} (mV)','V_{init} (mV)','V_{thresh} (mV)', 'Vslope_{max} (mV)'...
    ,'V_{half} (mV)','Spike_{amp} (mV)',... 'FaceColor','k'
    'AHP_{max} (mV)', 'Spike init (ms)','Spike_{rise} (ms)','Spike_{base width} (ms)','Spike_{half width} (ms)'};
for i=1:size(rk_all,1)
    subplot(2,9,i)
     histogram(rk_all(i,:),5,'FaceColor',[0.5 0.5 0.5]);box off;ylabel('Cells');
     hold on;
    title(str{i});ylim([0 60]);hold on;
    plot(nanmedian(rk_all(i,:)),55,'kv');
end

%% Plot spike diff waveform as imagesc
fig9=figure;set(fig9, 'Position', [200, 200, 800, 800]);set(gcf,'color','w');
for i=1:size(rk_all,2)
    subplot(8,8,i)
    try
     imagesc(trace_rk(st_rk(i):st_rk(i)+5*srF,i)');box off;
     hold on;
    catch
        imagesc(NaN);
    end
end
%% Plot diff spike waveform as imagesc
fig9=figure;set(fig9, 'Position', [200, 200, 800, 800]);set(gcf,'color','w');
for i=1:size(rk_all,2)
    subplot(8,8,i)
    try
     imagesc(diff(trace_rk(st_rk(i):st_rk(i)+5*srF,i))');box off;
     hold on;
    catch
        imagesc(NaN);
    end
end
%% Plot Frequcny 
plot_spike_freq(Ephys,rk);
%% Perform PCA and dip test
%remove NaNs
% reorder rk_all
rk_all=rk_all';
rk_all(find(sum(isnan(rk_all),2)>0),:)=[];
%% Run PCA
[coeff,score,latent,~,explained,mu] = pca([zscore(rk_all)]);
figure;bar(coeff(:,1),'k');set(gcf,'color','w');box off;ylabel('PC1 weight');xlabel('Feature number');set(gca,'FontSize',14);

%% Run dipt test
dip_test_SW(rk_all,1,{'PC1','PC2','PC3'})








%% Instrinisc properties Ntsr1 vs retro in the Ntsr mouse line
%retro
[rmp_retro maxsp_retro rheo_retro rin_retro tau_retro sag_retro trace_retro st_retro] = passive_readout(Ephys,retro_k);
%ntsr
[rmp_ntsr maxsp_ntsr rheo_ntsr rin_ntsr tau_ntsr sag_ntsr trace_ntsr st_ntsr] = passive_readout(Ephys,ntsr_k);
rmp_1= [rmp_retro(1:12)' rmp_ntsr'];
maxsp_1=[maxsp_retro(1:12)' maxsp_ntsr'];
rheo_1=[rheo_retro(1:12)' rheo_ntsr'];
rin_1=[rin_retro(1:12)' rin_ntsr'];
tau_1=[tau_retro(1:12)' tau_ntsr'];
sag_1=[sag_retro(1:12)' sag_ntsr'];
%% Reading out further active spike properties using Pandora
active_ntsr=sp_parameters_pandora(trace_ntsr,2);
active_retro=sp_parameters_pandora(trace_retro,2);
%% Plot exampe IVs of retro and NTSR1
cnr=11;step=[1:2:6 8:2:13];
ov_min=-150;ov_max=50;
temp=[];
temp=find(retro_k==1);
fig4=figure;set(fig4, 'Position', [200, 800, 400, 200]);set(gcf,'color','w');
subplot(1,2,1)
plot(Ephys(temp(cnr)).IV.traces(:,step),'Color','k','LineWidth',1);set(gca,'box','off');
hold on;ylim([ov_min-10 ov_max]);
hold on;text(-550*srF,Ephys(temp(cnr)).IV.RMP,[num2str(Ephys(temp(cnr)).IV.RMP),'mV'],'FontSize',9);
hold on;title('CPN');axis off;
subplot(1,2,2)
temp=[];
temp=find(ntsr_k==1);
plot(Ephys(temp(cnr)).IV.traces(:,step),'Color','r','LineWidth',1);set(gca,'box','off');
hold on;ylim([ov_min-10 ov_max]);
hold on;text(-550*srF,Ephys(temp(cnr)).IV.RMP,[num2str(Ephys(temp(cnr)).IV.RMP),'mV'],'FontSize',9,'Color','r');title('Ntsr1','Color','r');
axis off;
%Scale bar
 scale_x= 200;
 scale_y= 40;
 %scale barx
 hold on;x1= (1250*srF)-(100*srF);x2=1250*srF;p1=plot([x1 x2],[ov_min-10 ov_min-10],'-','Color','k','LineWidth',1.5);
 %scale bary
 hold on;y2= (ov_min-10)+scale_y;y1=(ov_min-10);p2=plot([x2 x2],[y1 y2],'-','Color','k','LineWidth',1.5); 
%% Test rmp, etc
close all;
%RMP
p1=paired_plot(rmp_1,0,{'k','r'});xticklabels({'CPN','NtsR1'});ylabel('RMP (mV)');yticks([-80:5:-60]);set(gca,'FontSize',10);
%Input resistance Rin
p2=paired_plot(rin_1,0,{'k','r'});xticklabels({'CPN','NtsR1'});ylabel('Input resistance (mOhm)');set(gca,'FontSize',10);
%Tau 
p3=paired_plot(tau_1,0,{'k','r'});xticklabels({'CPN','NtsR1'});ylabel('Tau (ms)');set(gca,'FontSize',10);
%Maximum spike rate
p4=paired_plot(maxsp_1,0,{'k','r'});xticklabels({'CPN','NtsR1'});ylabel('Max spike number');yticks([0:5:30]);set(gca,'FontSize',10);
%Rheobase
p5=paired_plot(rheo_1,1,{'k','r'});xticklabels({'CPN','NtsR1'});ylabel('Rheobase (pA)');yticks([0:50:200]);set(gca,'FontSize',10);
%Sag Ratio
p6=paired_plot(sag_1,1,{'k','r'});xticklabels({'CPN','NtsR1'});ylabel('Sag Ratio');set(gca,'FontSize',10);
%% Paired timing of first two spikes
p7=paired_plot(tim_1,1);xticklabels({'CPN','NtsR1'});ylabel('initial spikes');set(gca,'FontSize',10);
%% Plotting paired using the second spike always (if possible) 
close all;
p_val=paired_subplot(active_retro([1 2 3 4 5 6 7 8 10 11 14 15],:),active_ntsr([1 2 3 4 5 6 7 8 10 11 14 15],:),1);
%% plot imagsec from pvalues
p_all=horzcat(p1,p2,p3,p4,p5,p6,p_val)
fig3= figure;set(fig3, 'Name', '');set(fig3, 'Position', [200, 300, 600, 120]);set(gcf,'color','w');
str={'RMP','Rin','Tau','max nr spk','Rheo','Sag','V_{min}','V_{peak}','V_{init}','V_{thresh}', 'Vslope_{max}','V_{half}','Spike_{amplitude}',...
    'AHP_{max}', 'Spike init','Spike_{rise}','Spike_{base width}','Spike_{half width}'};
p_all(find(p_all>0.05))=1;
imagesc(p_all);
xticks([1:18])
xticklabels(str);
yticklabels({''});xtickangle(45);[cmap]=buildcmap('mw');colormap(cmap);colorbar;set(gca,'FontSize',10);
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
%% Plot exampe epsc / ipsc of retro and NTSR1
cnr=1;
ov_min=-150;ov_max=800;
temp=[];
temp=find(retro_cs==1);
fig4=figure;set(fig4, 'Position', [200, 800, 400, 200]);set(gcf,'color','w');
subplot(1,2,1)
plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,1),'Color','r','LineWidth',1);set(gca,'box','off');
hold on;plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,2),'Color','b','LineWidth',1);set(gca,'box','off');hold on;ylim([ov_min-10 ov_max]);
hold on;plot([0.25*sr 0.25*sr],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');

hold on;title('CPN');axis off;
subplot(1,2,2)
temp=[];
temp=find(ntsr_cs==1);
plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,1),'Color','r','LineWidth',1);set(gca,'box','off');
hold on;plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,2),'Color','b','LineWidth',1);set(gca,'box','off');
hold on;ylim([ov_min-10 ov_max]);title('Ntsr1','Color','k');
hold on;plot([0.25*sr 0.25*sr],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
axis off;
%Scale bar
 scale_x= 200;
 scale_y= 40;
 %scale barx
 hold on;x1= (1250*srF)-(100*srF);x2=1250*srF;p1=plot([x1 x2],[ov_min-10 ov_min-10],'-','Color','k','LineWidth',1.5);
 %scale bary
 hold on;y2= (ov_min-10)+scale_y;y1=(ov_min-10);p2=plot([x2 x2],[y1 y2],'-','Color','k','LineWidth',1.5);
 %% 
 cnr=4;
ov_min=-150;ov_max=30;
temp=[];
temp=find(retro_cs==1);
fig4=figure;set(fig4, 'Position', [200, 800, 400, 200]);set(gcf,'color','w');
subplot(1,2,1)
plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,1),'Color',[0 0.5 0.5],'LineWidth',1);set(gca,'box','off');
hold on;plot([0.25*sr 0.25*sr],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
hold on;title('CP');axis off;ylim([ov_min ov_max])

subplot(1,2,2)
temp=[];
temp=find(ntsr_cs==1);
plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,1),'Color',[0.7 0 0.4],'LineWidth',1);set(gca,'box','off');
title('CT','Color','k');
hold on;plot([0.25*sr 0.25*sr],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
axis off;ylim([ov_min ov_max])
%Scale bar
 scale_x= 200;
 scale_y= 40;
 %scale barx
 hold on;x1= (1250*srF)-(100*srF);x2=1250*srF;p1=plot([x1 x2],[ov_min-10 ov_min-10],'-','Color','k','LineWidth',1.5);
 %scale bary
 hold on;y2= (ov_min-10)+scale_y;y1=(ov_min-10);p2=plot([x2 x2],[y1 y2],'-','Color','k','LineWidth',1.5);
%% Plot EPSC and IPSC retro vs NtsR1
%EPSC
%color_id={[0.7 0 0.4],[0 0.5 0.5]};
paired_plot(epsc_1,0,{[0 0.5 0.5],[0.7 0 0.4]});xticklabels({'CP','CT'});ylabel('EPSC peak (pA)');set(gca,'FontSize',10);
%IPSC
paired_plot(ipsc_1,0,{[0 0.5 0.5],[0.7 0 0.4]});xticklabels({'CP','CT'});ylabel('IPSC peak (pA)');set(gca,'FontSize',10);
%E_I
paired_plot(e_i_1,0,{[0 0.5 0.5],[0.7 0 0.4]});xticklabels({'CP','CT'});ylabel('E / I Ratio');set(gca,'FontSize',10);
%% Read out epsp retro vs NtsR1
temp=[];
for i=1:length(find(retro_k==1));
    temp=find(retro_k==1);
    retro_epsp(i)=max(abs(Ephys(temp(i)).train_p(1,:)));   
end
temp=[];
for i=1:length(find(ntsr_k==1));
    temp=find(ntsr_k==1);
    ntsr_epsp(i)=max(abs(Ephys(temp(i)).train_p(1,:)));
end
%% Plot EPSP (in CC mode) for Ntsr1 vs CPN
epsp_1=[[retro_epsp(9:end)]' [ntsr_epsp]'];
paired_plot(epsp_1,0,{'k','r'});xticklabels({'CPN','NtsR1'});ylabel('EPSP peak (mV)');set(gca,'FontSize',10);
%% Plot exampe EPSP Ntsr1 and CPN
cnr=12;
ov_min=5;ov_max=20;
temp=[];
temp=find(retro_k==1);
fig4=figure;set(fig4, 'Position', [200, 300, 400, 200]);set(gcf,'color','w');
subplot(1,2,1)
plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,3),'Color','k','LineWidth',1);set(gca,'box','off');
hold on;plot([0.25*sr 0.25*sr],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
hold on;title('CPN');axis off;

cnr=8;
subplot(1,2,2)
temp=[];
temp=find(ntsr_k==1);
plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,2),'Color','r','LineWidth',1);set(gca,'box','off');
hold on;ylim([ov_min-10 ov_max]);title('Ntsr1','Color','k');
hold on;plot([0.25*sr 0.25*sr],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
axis off;
%Scale bar
 scale_x= 200;
 scale_y= 2.5;
 %scale barx
 hold on;x1= (1250*srF)-(100*srF);x2=1250*srF;p1=plot([x1 x2],[ov_min-10 ov_min-10],'-','Color','k','LineWidth',1.5); %scale bary
 hold on;y2= (ov_min-10)+scale_y;y1=(ov_min-10);p2=plot([x2 x2],[y1 y2],'-','Color','k','LineWidth',1.5); 
%% Time to peak ex and in for retro cells
%Ntsr1 mouse line, Cs-gluc, retro cells
rk_cs = cell_selecter(Ephys, 'label',1, 'sol',2);
temp=[];
temp=find(rk_cs==1);
close all
for i=1:length(temp)
[t_ex(i)]=time_to_peak(Ephys(temp(i)).sub_traces_train(:,1),[3000:30:5000],[5000:30:6000],20,0);
end
close all;
for i=1:length(temp)
[t_in(i)]=time_to_peak(Ephys(temp(i)).sub_traces_train(:,2),[3000:30:5000],[5000:30:6000],20,1);
end
close all;
%% Combine and remove cells where EX or IN is missing (NaN)
t_ein=[];
t_ein=[t_ex' t_in'];
[gh gm]=find(isnan(t_ein));
t_ein(unique(gh),:)=[];
%Plot paired and compare paired
paired_plot(t_ein,1,{'r','b'});xticklabels({'EX','IN'});ylabel('Time to peak (ms)');set(gca,'FontSize',10);
%% Time to peak ex and in for ntsr cell
temp=[];
temp=find(ntsr_cs==1);
close all
for i=1:length(temp)
[t_ex_ntsr(i)]=time_to_peak(Ephys(temp(i)).sub_traces_train(:,1),[3000:30:5000],[5000:30:6000],20,0);
end
for i=1:length(temp)
[t_in_ntsr(i)]=time_to_peak(Ephys(temp(i)).sub_traces_train(:,2),[3000:30:5000],[5000:30:6000],20,1);
end






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
epsc_2=[[retro_epsc_pv]' [pv_epsc]'];
ipsc_2=[retro_ipsc_pv' pv_ipsc'];
e_i_retro_pv=retro_epsc_pv./retro_ipsc_pv;
e_i_pv=pv_epsc./pv_ipsc;
e_i_2=[e_i_retro_pv' e_i_pv'];

%% 

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
%% Plot EPSC and IPSC retro vs PV
%EPSC
paired_plot(epsc_2,0,{[0 0.5 0.5],[0.8500, 0.3250, 0.0980]});xticklabels({'CP','PV'});ylabel('EPSC peak (pA)');set(gca,'FontSize',10);
%IPSC
paired_plot(ipsc_2,0,{'k','m'});xticklabels({'CPN','PV'});ylabel('IPSC peak (pA)');set(gca,'FontSize',10);
%E_I
paired_plot(e_i_2,0,{'k','m'});xticklabels({'CPN','PV'});ylabel('E / I Ratio');set(gca,'FontSize',10);

%% Plot all 5 iterations
temp=[];
for i=1:length(find(retro_cs_pv==1));
    temp=find(retro_cs_pv==1);
    retro_epsc_pv_5(:,i)=abs(Ephys(temp(i)).train_n(:,1));
    retro_ipsc_pv_5(:,i)=abs(Ephys(temp(i)).train_p(:,2));
end
temp=[];
for i=1:length(find(pv_cs==1));
    temp=find(pv_cs==1);
    pv_epsc_5(:,i)=abs(Ephys(temp(i)).train_n(:,1));
    pv_ipsc_5(:,i)=abs(Ephys(temp(i)).train_p(:,2));
end
figure;plot(retro_epsc_pv_5,'-ok');hold on;plot(pv_epsc_5,'-om')
%% Using all Interneurons with K
temp=[];
for i=1:length(find(ink==1));
    temp=find(ink==1);
    ink_epsc(i)=max(abs(Ephys(temp(i)).train_n(1,:)));    
end
temp=[];
for i=1:length(find(ink==1));
    temp=find(ink==1);
    ink_epsp(i)=max(abs(Ephys(temp(i)).train_p(1,:)));    
end
temp=[];
for i=1:length(find(rk==1));
    temp=find(rk==1);
    rk_epsc(i)=max(abs(Ephys(temp(i)).train_n(1,:)));    
end

temp=[];
for i=1:length(find(rk==1));
    temp=find(rk==1);
    rk_epsp(i)=max(abs(Ephys(temp(i)).train_p(1,:)));    
end
%% read out fast spiking interneurons
plot_intrinsic(Ephys,ink,7,srF,'m',1);
%% 
plot_intrinsic(Ephys,rk,7,srF,'k',1);

%% 
idx_fs=[1 3 5 6 9 10 11 12 13 14 16 17 18 20 27 28 29 31 32 34 35 36 37 38];
idx_nfs=[2 4 8 19 23 24 33];
idx_rks=[14:51];
fs_epsp=ink_epsp(idx_fs);
nfs_epsp=ink_epsp(idx_nfs);
rks_epsp=rk_epsp([idx_rks]);
%% Read out ephys IN vs CPN
[rmp_ink maxsp_ink rheo_ink rin_ink tau_ink sag_ink trace_ink st_ink] = passive_readout(Ephys,ink);
[rmp_rk maxsp_rk rheo_rk rin_rk tau_rk sag_rk trace_rk st_rk] = passive_readout(Ephys,rk);
%% Read out active intrinsic properties
active_ink=sp_parameters_pandora(trace_ink,2);
ink_all=vertcat(rmp_ink,rin_ink,tau_ink,maxsp_ink,rheo_ink,sag_ink,active_ink([1 2 3 4 5 6 7 8 10 11 14 15],:));
fig9=figure;set(fig9, 'Position', [200, 200, 1400, 400]);set(gcf,'color','w');
str={'RMP (mV)','R_{IN}(M\Omega)','Tau (ms)','Max Nr. Spikes','Rheobase (pA)','Sag ratio','V_{min} (mV)','V_{peak} (mV)','V_{init} (mV)','V_{thresh} (mV)', 'Vslope_{max} (mV)'...
    ,'V_{half} (mV)','Spike_{amp} (mV)',... 'FaceColor','k'
    'AHP_{max} (mV)', 'Spike init (ms)','Spike_{rise} (ms)','Spike_{base width} (ms)','Spike_{half width} (ms)'};
for i=1:size(ink_all,1)
    subplot(2,9,i)
     histogram(ink_all(i,:),5,'FaceColor',[0.5 0.5 0.5]);box off;ylabel('Cells');
     hold on;
    title(str{i});ylim([0 60]);hold on;
    plot(nanmedian(ink_all(i,:)),55,'kv');
end
%% Compare Rin between CPNs and all FS-IN
par1=[];
groups_idx=[];
%par1=vertcat(rin_rk(idx_rks)'./tau_rk(idx_rks)',rin_ink(idx_fs)'./tau_ink(idx_fs)')
%par1=vertcat(active_rk(3,idx_rks)',active_ink(3,idx_fs)')
par1=vertcat(active_rk(3,idx_rks)',active_ink(3,idx_fs)')
groups_idx=vertcat(ones(length(rin_rk(idx_rks)),1)*1,ones(length(rin_ink(idx_fs)),1)*2)
groups_idx(find(isnan(par1)))=[];
par1(find(isnan(par1)))=[];
[statsout] = barplot_sw(par1,groups_idx,{'','R_{IN}(M\Omega)'});xticklabels({'CPN','FS-IN'});xtickangle(45)
%% Rin vs EPSP for IN vs CPN
figure;scatter(fs_epsp,rin_ink(idx_fs)); hold on;scatter(rks_epsp,rin_rk(idx_rks));xlim([0 20])
%% Plot FS cell vs non FS GAD+
 %% Represantative example traces for nCPN CP and IN
cnr=16;
ov_min=5;ov_max=100;
temp=[];
temp=find(ink==1);
fig4=figure;set(fig4, 'Position', [200, 300, 400, 200]);set(gcf,'color','w');
subplot(1,2,1)
plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,2),'Color','m','LineWidth',1);set(gca,'box','off');
hold on;plot([0.25*sr 0.25*sr],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
hold on;title('FS-GAD+');axis off;
cnr=19;
subplot(1,2,2)
temp=[];
temp=find(ink==1);
plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,2),'Color','m','LineWidth',1);set(gca,'box','off');

hold on;ylim([ov_min-10 ov_max]);title('nonFS-GAD+','Color','k');
hold on;plot([0.25*sr 0.25*sr],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
axis off;
%Scale bar
 scale_x= 200;
 scale_y= 5;
 %scale barx
 hold on;x1= (1250*srF)-(100*srF);x2=1250*srF;p1=plot([x1 x2],[ov_min-10 ov_min-10],'-','Color','k','LineWidth',1.5);
 %scale bary
 hold on;y2= (ov_min-10)+scale_y;y1=(ov_min-10);p2=plot([x2 x2],[y1 y2],'-','Color','k','LineWidth',1.5); 

%% Compare EPSP FS-IN and non FS-IN
par1=[];groups_idx=[];
par1=vertcat(fs_epsp',nfs_epsp')
groups_idx=vertcat(ones(length(fs_epsp),1)*1,ones(length(nfs_epsp),1)*2)
groups_idx(find(isnan(par1)))=[];
par1(find(isnan(par1)))=[];
[statsout] = barplot_sw(par1,groups_idx,{'','EPSP peak (mV)'});xticklabels({'FS-GAD+','nonFS-GAD+'});xtickangle(45)

%% 




%% Plot CPN vs non labelled in  NTSR1 mice
temp=[];
for i=1:length(find(non_ntsr_lab==1));
    temp=find(non_ntsr_lab==1);
    retro_k_epsc_ntsr(i)=max(abs(Ephys(temp(i)).train_n(1,:)));    
end

 temp=[];
for i=1:length(find(non_ntsr==1));
    temp=find(non_ntsr==1);
    non_k_epsc_ntsr(i)=max(abs(Ephys(temp(i)).train_n(1,:)));    
end
 non_ntsr_retro_c=[retro_k_epsc_ntsr' non_k_epsc_ntsr'];
 %% 
 temp=[];
for i=1:length(find(non_ntsr_lab==1));
    temp=find(non_ntsr_lab==1);
    retro_k_epsp_ntsr(i)=max(abs(Ephys(temp(i)).train_p(1,:)));    
end

 temp=[];
for i=1:length(find(non_ntsr==1));
    temp=find(non_ntsr==1);
    non_k_epsp_ntsr(i)=max(abs(Ephys(temp(i)).train_p(1,:)));    
end
 non_ntsr_retro_p=[retro_k_epsp_ntsr' non_k_epsp_ntsr'];
 
 %% 
 paired_plot(non_ntsr_retro_c,0,{[0 0.5 0.5],[0.5 0.5 0.5]});xticklabels({'CP','CC'});ylabel('EPSC peak (pA)');set(gca,'FontSize',10);
 %% 
  paired_plot(non_ntsr_retro_p,0,{[0 0.5 0.5],[0.5 0.5 0.5]});xticklabels({'CP','CC'});ylabel('EPSP peak (mV)');set(gca,'FontSize',10);
 
  %% plot exmaple CP vs CC
  cnr=9;
ov_min=-150;ov_max=30;
temp=[];
temp=find(non_ntsr_lab==1);
fig4=figure;set(fig4, 'Position', [200, 800, 400, 200]);set(gcf,'color','w');
subplot(1,2,1)
plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,1),'Color',[0 0.5 0.5],'LineWidth',1);set(gca,'box','off');
hold on;plot([0.25*sr 0.25*sr],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
hold on;title('CP');
ylim([ov_min ov_max]);axis off

subplot(1,2,2)
temp=[];
temp=find(non_ntsr==1);
plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,1),'Color',[0.5 0.5 0.5],'LineWidth',1);set(gca,'box','off');
title('CC','Color','k');
hold on;plot([0.25*sr 0.25*sr],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
ylim([ov_min ov_max]);axis off
%Scale bar
 scale_x= 200;
 scale_y= 40;
 %scale barx
 hold on;x1= (1250*srF)-(100*srF);x2=1250*srF;p1=plot([x1 x2],[ov_min-10 ov_min-10],'-','Color','k','LineWidth',1.5);
 %scale bary
 hold on;y2= (ov_min-10)+scale_y;y1=(ov_min-10);p2=plot([x2 x2],[y1 y2],'-','Color','k','LineWidth',1.5); 
  
 %% 
 %% Plot IV for all non labelled cells
plot_intrinsic(Ephys,non_ntsr,4,srF,[0.5 0.5 0.5],1);
plot_intrinsic(Ephys,non_ntsr_lab,4,srF,[0.5 0.5 0.5],1);
 %% Instrinisc properties Ntsr1 vs retro in the Ntsr mouse line
%retro
[rmp_retro maxsp_retro rheo_retro rin_retro tau_retro sag_retro trace_retro st_retro] = passive_readout(Ephys,non_ntsr_lab);
%non lablled
[rmp_nl maxsp_nl rheo_nl rin_nl tau_nl sag_nl trace_nl st_nl] = passive_readout(Ephys,non_ntsr);


rmp_1= [rmp_retro' rmp_nl'];
maxsp_1=[maxsp_retro' maxsp_nl'];
rheo_1=[rheo_retro' rheo_nl'];
rin_1=[rin_retro' rin_nl'];
tau_1=[tau_retro' tau_nl'];
sag_1=[sag_retro' sag_nl'];
%% Reading out further active spike properties using Pandora
active_nl=sp_parameters_pandora(trace_nl,2);
active_retro=sp_parameters_pandora(trace_retro,2);
 %% Test rmp, etc
close all;
%RMP
p1=paired_plot(rmp_1,0,{'k','m'});xticklabels({'CPN','NtsR1'});ylabel('RMP (mV)');yticks([-80:5:-60]);set(gca,'FontSize',10);
%Input resistance Rin
p2=paired_plot(rin_1,0,{'k','m'});xticklabels({'CPN','NtsR1'});ylabel('Input resistance (mOhm)');set(gca,'FontSize',10);
%Tau 
p3=paired_plot(tau_1,0,{'k','m'});xticklabels({'CPN','NtsR1'});ylabel('Tau (ms)');set(gca,'FontSize',10);
%Maximum spike rate
p4=paired_plot(maxsp_1,0,{'k','m'});xticklabels({'CPN','NtsR1'});ylabel('Max spike number');yticks([0:5:30]);set(gca,'FontSize',10);
%Rheobase
p5=paired_plot(rheo_1,0,{'k','m'});xticklabels({'CPN','NtsR1'});ylabel('Rheobase (pA)');yticks([0:50:200]);set(gca,'FontSize',10);
%Sag Ratio
p6=paired_plot(sag_1,0,{'k','m'});xticklabels({'CPN','NtsR1'});ylabel('Sag Ratio');set(gca,'FontSize',10);
%% Paired timing of first two spikes
p7=paired_plot(tim_1,1);xticklabels({'CPN','NtsR1'});ylabel('initial spikes');set(gca,'FontSize',10);
%% Plotting paired using the second spike always (if possible) 
close all;
p_val=paired_subplot(active_retro([1 2 3 4 5 6 7 8 10 11 14 15],1:8),active_nl([1 2 3 4 5 6 7 8 10 11 14 15],:),1);


%% plot imagsec from pvalues
p_all=horzcat(p1,p2,p3,p4,p5,p6,p_val)
fig3= figure;set(fig3, 'Name', '');set(fig3, 'Position', [200, 300, 600, 120]);set(gcf,'color','w');
str={'RMP','Rin','Tau','max nr spk','Rheo','Sag','V_{min}','V_{peak}','V_{init}','V_{thresh}', 'Vslope_{max}','V_{half}','Spike_{amplitude}',...
    'AHP_{max}', 'Spike init','Spike_{rise}','Spike_{base width}','Spike_{half width}'};
p_all(find(p_all>0.05))=1;
imagesc(p_all);
xticks([1:18])
xticklabels(str);
yticklabels({''});xtickangle(45);[cmap]=buildcmap('mw');colormap(cmap);colorbar;set(gca,'FontSize',10);
%% 


%% TTX experiments
ttx_id=[206 207 208 209 214];
%% Plot before and after TTX
cnr=ttx_id(3);
ov_min=-20;ov_max=600;
fig4=figure;set(fig4, 'Position', [200, 800, 400, 200]);set(gcf,'color','w');
subplot(1,2,1)
plot(Ephys(cnr).sub_traces_train(1:1*sr,2),'Color',[0.5 0.5 0.5],'LineWidth',1.2);set(gca,'box','off');
hold on;plot(Ephys(cnr).sub_traces_train(1:1*sr,3),'Color','r','LineWidth',1.2);set(gca,'box','off');
axis off;
hold on;plot([0.25*sr 0.25*sr],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
%scale barx
 hold on;x1= (1250*srF)-(100*srF);x2=1250*srF;p1=plot([x1 x2],[ov_min-10 ov_min-10],'-','Color','k','LineWidth',1.5);
 %scale bary
 hold on;y2= (ov_min-10)+scale_y;y1=(ov_min-10);p2=plot([x2 x2],[y1 y2],'-','Color','k','LineWidth',1.5);
subplot(1,2,2)
cnr=ttx_id(3);
ov_min=-100;ov_max=20;
plot(Ephys(cnr).sub_traces_train(1:1*sr,1),'Color',[0.5 0.5 0.5],'LineWidth',1.2);set(gca,'box','off');
hold on;plot(Ephys(cnr).sub_traces_train(1:1*sr,4),'Color','r','LineWidth',1.2);set(gca,'box','off');
axis off;
hold on;plot([0.25*sr 0.25*sr],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
%Scale bar
 scale_x= 200;
 scale_y= 40;
 %scale barx
 hold on;x1= (1250*srF)-(100*srF);x2=1250*srF;p1=plot([x1 x2],[ov_min-10 ov_min-10],'-','Color','k','LineWidth',1.5);
 %scale bary
 hold on;y2= (ov_min-10)+scale_y;y1=(ov_min-10);p2=plot([x2 x2],[y1 y2],'-','Color','k','LineWidth',1.5);
%% Quantification for IPSC
 temp=[];ttx_p=[];
 temp=ttx_id([1 3 4]);
for i=1:length(ttx_id([1 3 4]));
    ttx_p(i,:)=[max(abs(Ephys(temp(i)).train_p(:,2))) max(abs(Ephys(temp(i)).train_p(:,3)))];    
end
ttx_214=[max(abs(Ephys(214).train_p(:,1))) max(abs(Ephys(214).train_p(:,2)))];
ttx_all_p=[ttx_p; ttx_214]
ttx_all_p./max(ttx_all_p,[],2);
p1=paired_plot(ttx_all_p,1,{'[0.5 0.5 0.5]','r'});xticklabels({'Before','After'});ylabel('Peak IPSC (pA)');set(gca,'FontSize',10);
%% Quantification for EPSC
 temp=[];ttx_n=[];
 temp=ttx_id([1 2 3 4]);
for i=1:length(temp)
    ttx_n(i,:)=[max(abs(Ephys(temp(i)).train_n(:,1))) max(abs(Ephys(temp(i)).train_n(:,4)))];    
end
ttx_214=[max(abs(Ephys(214).high_n(:,1))) max(abs(Ephys(214).high_n(:,4)))];
ttx_all_n=[ttx_n; ttx_214]

p1=paired_plot(ttx_all_n,1,{'[0.5 0.5 0.5]','r'});xticklabels({'Before','After'});ylabel('Peak EPSC (pA)');set(gca,'FontSize',10);
%% Compare PV EPSP (spike no spike) with CPs
temp=[];
for i=1:length(find(retro_k_pv==1));
    temp=find(retro_k_pv==1);
    r_epsp_pv(i)=max(abs(max(Ephys(temp(i)).train_p(:)))); 
     try
   r_epsp_pvh(i)=max(abs(max(Ephys(temp(i)).high_p(:))));  
    catch
    r_epsp_pvh(i)=NaN;
    end
end
temp=[];
 temp=find(pv_k==1);
 temp(find([Ephys(pv_k).wc]==0))=[]
for i=1:length(temp);
    pv_epsp(i)=max(abs(max(Ephys(temp(i)).train_p(:))));
    try
    pv_epsp_h(i)=max(abs(max(Ephys(temp(i)).high_p(:))));  
    catch
     pv_epsp_h(i)=NaN;
    end
end
retro_epsp_pv=nanmax(cat(1,r_epsp_pv,r_epsp_pvh));
pv_epsp_pv=nanmax(cat(1,pv_epsp,pv_epsp_h));
%% Compare EPSP between CP and PV
par=[retro_epsp_pv pv_epsp_pv];
g1=1:length(retro_epsp_pv);
g2=length(retro_epsp_pv)+1:length(pv_epsp_pv)+length(retro_epsp_pv);
[statsout]=dual_barplot(par,g1,g2,2); ylabel('EPSP (mV)');xticks([1:1:2]);xticklabels({'CP','PV'});set(gca,'FontSize',10);
%% 
%% Plot exampe EPSP PV and CP


fig4=figure;set(fig4, 'Position', [200, 300, 400, 200]);set(gcf,'color','w');
ov_min=-20;ov_max=80;
subplot(1,2,1)
cnr=1;
temp=[];
temp=find(retro_k_pv==1);
plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,2),'Color',[0.5 0.5 0.5],'LineWidth',1);set(gca,'box','off');
hold on;ylim([ov_min-15 ov_max]);title('CP','Color','k');
hold on;plot([0.25*sr 0.25*sr],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
axis off;ylim([ov_min ov_max])

subplot(1,2,2)
cnr=4;

temp=[];
temp=find(pv_k==1);
plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,2),'Color',[0.8500, 0.3250, 0.0980],'LineWidth',1);set(gca,'box','off');
hold on;plot([0.25*sr 0.25*sr],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
hold on;title('PV');axis off;ylim([ov_min ov_max])



%Scale bar
 scale_x= 200;
 scale_y= 2.5;
 %scale barx
 hold on;x1= (1250*srF)-(100*srF);x2=1250*srF;p1=plot([x1 x2],[ov_min ov_min],'-','Color','k','LineWidth',1.5); %scale bary
 hold on;y2= (ov_min)+scale_y;y1=(ov_min);p2=plot([x2 x2],[y1 y2],'-','Color','k','LineWidth',1.5); 

 %% pie chart
