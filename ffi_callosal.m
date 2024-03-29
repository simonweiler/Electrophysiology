%% FEEDFORWARD INHIBITION
%which cell types spike upon callosal activation in the slice? : CPN, CT,
%CCi (Penk and unlablled CCs), PV, SST, VIP, GAD
%% 

save_folder='C:\Users\simonw\S1-V1 interaction Dropbox\Simon Weiler\Callosal_L6\SfN2023\Figures';
%% Load data structure 
str_L6    = 'D:\Postdoc_Margrie\Projects\Callosal\output';
folder_list = uipickfiles('FilterSpec',str_L6);
load(char(folder_list));
%sampling rate
srF=20;sr=20000;
%% EPSC ISPC timing comparison CP and CT
%readout all cells type for timw windowing epsc ipsc
  cp_psc1=cell_selecter(Ephys,'drugs',0,'label',1,'sol',2,'qualityinput',1);
  cp_psc2=cell_selecter(Ephys,'drugs',1,'label',1,'sol',2,'qualityinput',1);
  cp_psc=[];cp_psc=cp_psc1+cp_psc2;
%readout amplitude if EPSC and IPSC
[cp_amp] = readout_amp_update(Ephys,cp_psc ,1,2,1,2);
% Time to peak ex and in for retro cells using first pulse of long train 
%CPN
%% 
rsp_window=[];
rsp_window=[5001:1:6000];
fc_1=7;
[time_cp1] = timing_readout(Ephys, cp_psc, cp_amp(:,3), fc_1, 1,1,rsp_window);
%% 
[time_cp2] = timing_readout(Ephys, cp_psc, cp_amp(:,3), fc_1, 1,2,rsp_window);
%% 
[time_cp3] = timing_readout(Ephys, cp_psc, cp_amp(:,3), fc_1, 1,3,rsp_window);
%%  Exclude cells with 25 ms and take the average of the first pulse across 1 , 5 and 10 Hz CP
time_cp1(find(time_cp1(:,1)>25),1)=NaN;
time_cp2(find(time_cp2(:,1)>25),1)=NaN;
time_cp3(find(time_cp3(:,1)>25),1)=NaN;
ex_temp_cp=[time_cp1(:,1) time_cp2(:,1) time_cp3(:,1)];
time_cp=[nanmean([time_cp1(:,1) time_cp2(:,1) time_cp3(:,1)],2) nanmean([time_cp1(:,2) time_cp2(:,2) time_cp3(:,2)],2) (time_cp1(:,3))];
time_cpstd=[nanstd([time_cp1(:,1) time_cp2(:,1) time_cp3(:,1)],[],2) nanstd([time_cp1(:,2) time_cp2(:,2) time_cp3(:,2)],[],2) (time_cp1(:,3))];
cp_avg=[];cp_avg=time_cp(find(sum(~isnan(ex_temp_cp),2)>1),[1 2 3]);

%% Plot bar plot comparison
%only include cells that have more than 1 ms EPSC (sone have indee shorter
%ones which are then direct?)
data=[];data=cp_avg(find(cp_avg(:,1)>1),[1 2]);

cl={'r','b'};
paired_plot_box(data,cl);
xticklabels({'EX','IN'});ylabel('Onset Latency (ms)');
%statistics 
%test for normality: 
set(gca,'FontSize',11);set(gca,'TickDir','out');
%pvalue of paired signrank test: 
[p1]=signrank(time_cp(:,1),time_cp(:,2))
%% add asteriks
title('');text(1,20,'***','FontSize',18,'FontWeight','normal');
%% SAVE
cd(save_folder);saveas(gcf, 'CP_EPSC_IPSC_onset.pdf');

%% %readout all cells type for timw windowing epsc ipsc
  ct_psc1=cell_selecter(Ephys,'drugs',0,'label',4,'geno',6,'sol',2,'qualityinput',1);
  ct_psc2=cell_selecter(Ephys,'drugs',1,'label',4,'geno',6,'sol',2,'qualityinput',1);
  ct_psc=[];ct_psc=ct_psc1+ct_psc2;
%readout amplitude if EPSC and IPSC
[ct_amp] = readout_amp_update(Ephys,ct_psc ,1,2,1,2);
% Time to peak ex and in for retro cells using first pulse of long train 
%CPN
%% 
rsp_window=[];
rsp_window=[5001:1:6000];
%rsp_window=[25001:1:26000];

fc_1=7;
[time_ct1]=[];
[time_ct1] = timing_readout(Ephys, ct_psc, ct_amp(:,3), fc_1, 1,1,rsp_window);
%% 
[time_ct2] = timing_readout(Ephys, ct_psc, ct_amp(:,3), fc_1, 1,2,rsp_window);
%% 

[time_ct3] = timing_readout(Ephys, ct_psc, ct_amp(:,3), fc_1, 1,3,rsp_window);
%% Exclude cells with 25 ms and take the average of the first pulse across 1 , 5 and 10 Hz CT
time_ct1(find(time_ct1(:,1)>25),1)=NaN;
time_ct2(find(time_ct2(:,1)>25),1)=NaN;
time_ct3(find(time_ct3(:,1)>25),1)=NaN;
ex_temp_ct=[time_ct1(:,1) time_ct2(:,1) time_ct3(:,1)];
time_ct=[nanmean([time_ct1(:,1) time_ct2(:,1) time_ct3(:,1)],2) nanmean([time_ct1(:,2) time_ct2(:,2) time_ct3(:,2)],2) (time_ct1(:,3))];
time_ctstd=[nanstd([time_ct1(:,1) time_ct2(:,1) time_ct3(:,1)],[],2) nanstd([time_ct1(:,2) time_ct2(:,2) time_ct3(:,2)],[],2) (time_ct1(:,3))];
ct_avg=[];ct_avg=time_ct(find(sum(~isnan(ex_temp_ct),2)>1),[1 2 3]);

%% Plot CT EPSC IPSC onset delay
data=[];data=ct_avg(find(ct_avg(:,1)>1),[1 2]);
cl={'r','b'};
paired_plot_box(data,cl);
xticklabels({'EX','IN'});ylabel('Onset Latency (ms)');
%statistics 
%test for normality: 
set(gca,'FontSize',11);set(gca,'TickDir','out');
%pvalue of paired signrank test: 
[p1]=signrank(time_ct(:,1),time_ct(:,2));ylim([0 20])
%% add asteriks
title('');text(1,20,'**','FontSize',18,'FontWeight','normal');
%% SAVE CT
cd(save_folder);saveas(gcf, 'CT_EPSC_IPSC_onset.pdf');



%% Read out all EPSPs across cell types
%Interneurons 
stim_type=1;
%PV
pv_k_cells=cell_selecter(Ephys,'drugs',0,'label',3,'sol',1,'qualityinput',1);
epsp_pv=[];
[epsp_pv] = readout_amp_epsp(Ephys,pv_k_cells,stim_type,sr);
disp([num2str(length(find(epsp_pv>45))) ' PV cells of ' num2str(length(epsp_pv)) ' spiked'])
%
% SOM
som_k_cells=cell_selecter(Ephys,'drugs',0,'label',6,'geno',8,'sol',1,'qualityinput',1);
[epsp_som]=[];
[epsp_som] = readout_amp_epsp(Ephys,som_k_cells,stim_type,sr);
disp([num2str(length(find(epsp_som>45))) ' SST cells of ' num2str(length(epsp_som)) ' spiked'])
%
% read out all GAD
in_k_cells=cell_selecter(Ephys,'drugs',0,'label',2,'sol',1,'qualityinput',1);
[epsp_gad]=[];
[epsp_gad] = readout_amp_epsp(Ephys,in_k_cells,stim_type,sr);
disp([num2str(length(find(epsp_gad>45))) ' GAD cells of ' num2str(length(epsp_gad)) ' spiked'])
%
% VIP
vip_k_cells=cell_selecter(Ephys,'drugs',0,'label',8,'geno',9,'sol',1,'optovariant',1);
[epsp_vip]=[];
[epsp_vip] = readout_amp_epsp(Ephys,vip_k_cells,stim_type,sr);
disp([num2str(length(find(epsp_vip>45))) ' VIP cells of ' num2str(length(epsp_vip)) ' spiked'])
%% Excitatry PNS
% CPN
cpn_k_cells=cell_selecter(Ephys,'drugs',0,'label',1,'sol',1,'qualityinput',1);
[epsp_cpn]=[];
[epsp_cpn] = readout_amp_epsp(Ephys,cpn_k_cells,stim_type,sr);
disp([num2str(length(find(epsp_cpn>45))) ' CP cells of ' num2str(length(epsp_cpn)) ' spiked'])

% Ntsr1
ntsr_k_cells=cell_selecter(Ephys,'drugs',0,'label',4,'geno',6,'sol',1,'qualityinput',1);
[epsp_ntsr]=[];
[epsp_ntsr] = readout_amp_epsp(Ephys,ntsr_k_cells,stim_type,sr);
disp([num2str(length(find(epsp_ntsr>45))) ' CT cells of ' num2str(length(epsp_ntsr)) ' spiked'])

% Penk
penk_k_cells=cell_selecter(Ephys,'drugs',0,'label',5,'geno',7,'sol',1,'qualityinput',1);
[epsp_penk]=[];
[epsp_penk] = readout_amp_epsp(Ephys,penk_k_cells,stim_type,sr); 
disp([num2str(length(find(epsp_penk>45))) ' Penk cells of ' num2str(length(epsp_penk)) ' spiked'])

% Unlabelled excitatry cell in ntsr1 animal (additional to Penk CCi) only K solution
cci_k_cells=cell_selecter(Ephys,'drugs',0,'label',0,'geno',6,'sol',1,'qualityinput',1);
[epsp_cci]=[];
[epsp_cci] = readout_amp_epsp(Ephys,cci_k_cells,stim_type,sr); 
disp([num2str(length(find(epsp_cci>45))) ' CCI cells of ' num2str(length(epsp_cci)) ' spiked'])
%% Display fraction per cell type
fig6= figure;set(fig6, 'Name', 'compare fraction spiking');set(fig6, 'Position', [200, 300, 300, 250]);set(gcf,'color','w');
gr_m=[(length(find(epsp_cpn>45))/length(epsp_cpn))*100  length(find(epsp_ntsr>45))/length(epsp_ntsr)*100 (length(find(epsp_penk>45))+length(find(epsp_cci>45)))/(length(epsp_penk)+length(epsp_cci))*100 ...
   NaN length(find(epsp_gad>45))/length(epsp_gad)*100 length(find(epsp_pv>45))/length(epsp_pv)*100 ...
    length(find(epsp_som>45))/length(epsp_som)*100 length(find(epsp_vip>45))/length(epsp_vip)*100];  
col_m=[([0 166 156]./256);([235 0 139]./256);([190 30 45]./256);([180 180 180]./256);([180 180 180]./256);([246 146 30]./256);([27 117 187]./256);([0 165 81]./256)];
for i=1:3
b=bar(i,gr_m(i));b.EdgeColor=col_m(i,:);b.FaceAlpha=0;b.LineWidth=1.5;
hold on;
text(i,round(gr_m(i),0),[num2str(round(gr_m(i),0)) '%'],'vert','bottom','horiz','center'); 
hold on;
end
hold on
for i=5:8
b=bar(i,gr_m(i));b.EdgeColor=col_m(i,:);b.FaceAlpha=0;b.LineWidth=1.5;
hold on;
text(i,round(gr_m(i),0),[num2str(round(gr_m(i),0)) '%'],'vert','bottom','horiz','center'); 
hold on;
end
xlim([0 9]);xticks([1:1:8]);ylabel('Percentage of cells firing APs');xticklabels({'CP','CT','CCi','','GAD','PV','SST','VIP'});xtickangle(45);set(gca,'FontSize',11);
ylim([0 50]);box off;h = gca;h.YAxis.Visible = 'off';title('Percentage of cells firing APs','FontSize',12,'FontWeight','normal');set(gca,'TickDir','out');
hold on;
text(1.5,40,'PNs','FontSize',12);text(5.75,40,'INs','FontSize',12);
%% SAVE
cd(save_folder);saveas(gcf, 'Fraction_spiking.pdf');





%% Read out max spike frequency for all cell types
[cpn_spk] = spike_counts(Ephys,cpn_k_cells);
[ntsr_spk] = spike_counts(Ephys,ntsr_k_cells);
[penk_spk] = spike_counts(Ephys,penk_k_cells);
[cci_spk] = spike_counts(Ephys,cci_k_cells);
[pv_spk] = spike_counts(Ephys,pv_k_cells);
[som_spk] = spike_counts(Ephys,som_k_cells);
[vip_spk] = spike_counts(Ephys,vip_k_cells);
[gad_spk] = spike_counts(Ephys,in_k_cells);
%% find threshold for fast spiker
threshold_FS=round((prctile([som_spk],25) + prctile([pv_spk],5) + prctile([cpn_spk ntsr_spk penk_spk cci_spk],95))/3)
%% scatter plot (spike freq vs EPSP amplitude)
%NOTE THAT ONE ST cell that fired has no max spiking hence its one less
%cell for SST

pointsize=100;
fig4=figure;set(fig4, 'Position', [200, 400, 550, 500]);set(gcf,'color','w');
% s1=scatter(cpn_spk,epsp_cpn,pointsize,'c','filled','MarkerFaceColor',([0 166 156]./256),'MarkerEdgeColor','w');hold on;
% s2=scatter(ntsr_spk,epsp_ntsr,pointsize,'m','filled','MarkerFaceColor',([235 0 139]./256),'MarkerEdgeColor','w');hold on;
% s3=scatter(penk_spk,epsp_penk,pointsize,'r','filled','MarkerFaceColor',([190 30 45]./256),'MarkerEdgeColor','w');hold on;
%s4=scatter(cci_spk,epsp_cci,pointsize,'r','filled','MarkerFaceColor',([190 30 45]./256),'MarkerEdgeColor','w');hold on;
s5=scatter(pv_spk,epsp_pv,pointsize,'y','filled','MarkerFaceColor',([246 146 30]./256),'MarkerEdgeColor','w');hold on;
s6=scatter(som_spk,epsp_som,pointsize,'b','filled','MarkerFaceColor',([27 117 187]./256),'MarkerEdgeColor','w');hold on;
s7=scatter(vip_spk,epsp_vip,pointsize,'g','filled','MarkerFaceColor',([0 165 81]./256),'MarkerEdgeColor','w');hold on;
s8=scatter(gad_spk,epsp_gad,pointsize,'k','filled','MarkerFaceColor',([0.5 0.5 0.5]),'MarkerEdgeColor','w');hold on;
xlim([-3 120]);ylim([-3 120]);
hold on;line([threshold_FS threshold_FS],[ylim],'linewidth',0.5,'Color','k','LineStyle','--');
hold on;line([xlim],[45 45],'linewidth',0.5,'Color','k','LineStyle','--');
% legend({'CPN','CT','CCi','CCi','PV','SST','VIP','GAD'},'Location','northwest');
legend({'PV','SST','VIP','GAD'},'Location','northwest');
legend boxoff
box off; ylabel('Light evoked response (mV)');xlabel('Max. spike frequency (Hz)');
text(threshold_FS+1,45+4,'spike threshold','FontSize',12);
text(10,120,'RS','FontSize',12);text(70,120,'FS','FontSize',12);
set(gca,'FontSize',12);set(gca,'TickDir','out');
%% SAVE
cd(save_folder);saveas(gcf, 'Spikes_vs_APFfreq_INs.pdf');



%% 
[rmp_cpn maxsp_cpn rheo_cpn rin_cpn tau_cpn sag_cpn trace_cpn spike_time_cpn] = passive_readout(Ephys,cpn_k_cells);
[rmp_cpn maxsp_ct rheo_ct rin_ct tau_ct sag_ct trace_ct spike_time_ct] = passive_readout(Ephys,ntsr_k_cells);
[rmp_penk maxsp_penk rheo_penk rin_penk tau_penk sag_penk trace_penk spike_time_penk] = passive_readout(Ephys,penk_k_cells);
[rmp_cci maxsp_cci rheo_cci rin_cci tau_cci sag_cci trace_cci spike_time_cci] = passive_readout(Ephys,penk_k_cells);
[rmp_pv maxsp_pv rheo_pv rin_pv tau_pv sag_pv trace_pv spike_time_pv] = passive_readout(Ephys,pv_k_cells);
[rmp_som maxsp_som rheo_som rin_som tau_som sag_som trace_som spike_time_som] = passive_readout(Ephys,som_k_cells);
[rmp_vip maxsp_vip rheo_vip rin_vip tau_vip sag_vip trace_vip spike_time_vip] = passive_readout(Ephys,vip_k_cells);
[rmp_gad maxsp_gad rheo_gad rin_gad tau_gad sag_gad trace_gad spike_time_gad] = passive_readout(Ephys,in_k_cells);
%% 
spike_nr=1;
active_pv=sp_parameters_pandora(trace_pv,spike_nr);
active_som=sp_parameters_pandora(trace_som,spike_nr);
active_vip=sp_parameters_pandora(trace_vip,spike_nr);
active_gad=sp_parameters_pandora(trace_gad,spike_nr);

%% check ephys propoerties for cell that spiked upon light activtaion and the non spiking cells
spike_thr=50;
max_sp_thr=30;

p1=ranksum([rmp_pv(find(epsp_pv>spike_thr)) rmp_som(find(epsp_som>spike_thr)) rmp_gad(find(epsp_gad>spike_thr))],...
    [rmp_som(find(som_spk<max_sp_thr)) rmp_gad(find(gad_spk<max_sp_thr)) rmp_vip(find(vip_spk<max_sp_thr))])
p2=ranksum([rin_pv(find(epsp_pv>spike_thr)) rin_som(find(epsp_som>spike_thr)) rin_gad(find(epsp_gad>spike_thr))],...
    [rin_som(find(som_spk<max_sp_thr)) rin_gad(find(gad_spk<max_sp_thr)) rin_vip(find(vip_spk<max_sp_thr))])
p3=ranksum([rheo_pv(find(epsp_pv>spike_thr)) rheo_som(find(epsp_som>spike_thr)) rheo_gad(find(epsp_gad>spike_thr))],...
    [rheo_som(find(som_spk<max_sp_thr)) rheo_gad(find(gad_spk<max_sp_thr)) rheo_vip(find(vip_spk<max_sp_thr))])
p4=ranksum([tau_pv(find(epsp_pv>spike_thr)) tau_som(find(epsp_som>spike_thr)) tau_gad(find(epsp_gad>spike_thr))],...
    [tau_som(find(som_spk<max_sp_thr)) tau_gad(find(gad_spk<max_sp_thr)) tau_vip(find(vip_spk<max_sp_thr))])
%% 
p_active=[];
for nr=1:18

p_active(nr)=ranksum([active_pv(nr,find(epsp_pv>spike_thr)) active_som(nr,find(epsp_som>spike_thr)) active_gad(nr,find(epsp_gad>spike_thr),1)],...
    [active_som(nr,find(som_spk<max_sp_thr)) active_gad(nr,find(gad_spk<max_sp_thr)) active_vip(nr,find(vip_spk<max_sp_thr))])
end
%% 

p_active=[];
for nr=1:18

p_active(nr)=ranksum([active_pv(nr,find(pv_spk>max_sp_thr)) active_som(nr,find(som_spk>max_sp_thr)) active_gad(nr,find(gad_spk>max_sp_thr),1)],...
    [active_som(nr,find(som_spk<max_sp_thr)) active_gad(nr,find(gad_spk<max_sp_thr)) active_vip(nr,find(vip_spk<max_sp_thr))])
end
%% 
p_active=[];
for nr=1:18

p_active(nr)=ranksum([active_pv(nr,find(epsp_pv>spike_thr)) active_som(nr,find(epsp_som>spike_thr)) active_gad(nr,find(epsp_gad>spike_thr),1)],...
    [active_pv(nr,find(epsp_pv<spike_thr)) active_som(nr,find(epsp_som<spike_thr)) active_gad(nr,find(epsp_gad<spike_thr),1)])
end

%% 

g1=[];g2=[];
p1=[];p2=[];
%combine ex and in
p1=[rin_pv(find(epsp_pv>spike_thr)) rin_som(find(epsp_som>spike_thr)) rin_gad(find(epsp_gad>spike_thr))];p2=[rin_som(find(som_spk<max_sp_thr)) rin_gad(find(gad_spk<max_sp_thr)) rin_vip(find(vip_spk<max_sp_thr))];
par=[];par=[p1 p2]';
g1=[];g1=ones(1,length(p1)); 
g2=[];g2=ones(1,length(p2))*2;
gro=[];gro=[g1 g2]';
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 300, 300]);set(gcf,'color','w');
%hold on;line([0 3],[0 0],'LineStyle',':','Color','k','LineWidth',1)
violins = violinplot(par, gro,'ViolinColor',[[0.3 0.3 0.3]; [0.6 0.6 0.6]]);box off;
xlim([0 3]);ylabel('Input Resistance');set(gca,'FontSize',12);
% h = gca;h.XAxis.Visible = 'off';
% set(gca,'xtick',[]);
xticklabels({'L6 FS INs','L6 NFS INs'});xtickangle(45);
%ylim([-1.2 1.2]);
%% 

g1=[];g2=[];
p1=[];p2=[];
%combine ex and in
p1=[tau_pv(find(epsp_pv>spike_thr)) tau_som(find(epsp_som>spike_thr)) tau_gad(find(epsp_gad>spike_thr))];p2=[tau_som(find(som_spk<max_sp_thr)) tau_gad(find(gad_spk<max_sp_thr)) tau_vip(find(vip_spk<max_sp_thr))];
par=[];par=[p1 p2]';
g1=[];g1=ones(1,length(p1)); 
g2=[];g2=ones(1,length(p2))*2;
gro=[];gro=[g1 g2]';
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 300, 300]);set(gcf,'color','w');
%hold on;line([0 3],[0 0],'LineStyle',':','Color','k','LineWidth',1)
violins = violinplot(par, gro,'ViolinColor',[[0.3 0.3 0.3]; [0.6 0.6 0.6]]);box off;
xlim([0 3]);ylabel('Tau');set(gca,'FontSize',12);
% h = gca;h.XAxis.Visible = 'off';
% set(gca,'xtick',[]);
xticklabels({'L6 FS INs','L6 NFS INs'});xtickangle(45);
%ylim([-1.2 1.2]);
%% 
nr=4
g1=[];g2=[];
p1=[];p2=[];
%combine ex and in
p1=[active_pv(nr,find(epsp_pv>spike_thr)) active_som(nr,find(epsp_som>spike_thr)) active_gad(nr,find(epsp_gad>spike_thr),1)];p2=[active_som(nr,find(som_spk<max_sp_thr)) active_gad(nr,find(gad_spk<max_sp_thr)) active_vip(nr,find(vip_spk<max_sp_thr))];
par=[];par=[p1 p2]';
g1=[];g1=ones(1,length(p1)); 
g2=[];g2=ones(1,length(p2))*2;
gro=[];gro=[g1 g2]';
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 300, 300]);set(gcf,'color','w');
%hold on;line([0 3],[0 0],'LineStyle',':','Color','k','LineWidth',1)
violins = violinplot(par, gro,'ViolinColor',[[0.3 0.3 0.3]; [0.6 0.6 0.6]]);box off;
xlim([0 3]);ylabel('Tau');set(gca,'FontSize',12);
% h = gca;h.XAxis.Visible = 'off';
% set(gca,'xtick',[]);
xticklabels({'L6 FS INs','L6 NFS INs'});xtickangle(45);
%ylim([-1.2 1.2]);

%% example injected current vs spike frequency

cnr=22;cnr2=26;
temp=[];temp=find(som_k_cells==1);stimvec_cpn=[];maxcount_spk=[];stimvec_cpn2=[];maxcount_spk2=[];
stimvec_cpn=Ephys(temp(cnr)).IV.stimvec;
maxcount_spk=Ephys(temp(cnr)).IV.spikecount;
stimvec_cpn2=Ephys(temp(cnr2)).IV.stimvec;
maxcount_spk2=Ephys(temp(cnr2)).IV.spikecount;

fig4=figure;set(fig4, 'Position', [200, 600, 300, 200]);set(gcf,'color','w');
     p1=plot(stimvec_cpn,maxcount_spk,'-o');p1.Color=([27 30 187]./256);p1.Color(4)=3/8; p1.MarkerSize=3;p1.MarkerFaceColor=([27 30 187]./256);p1.MarkerSize=5;
     hold on;set(gca,'box','off');hold on;
     p2=plot(stimvec_cpn2,maxcount_spk2,'-o');p2.Color=([27 117 187]./256);p2.Color(4)=3/8; p2.MarkerSize=3;p2.MarkerFaceColor=([27 117 187]./256);p2.MarkerSize=5;
     hold on;set(gca,'box','off'); 
ylabel('Spike frequency (Hz)');xlabel('Injected current (pA)')
xlim([0 500]);set(gca,'FontSize',12);
%% 

%show max spike frequcny trace
fig4=figure;set(fig4, 'Position', [200, 400, 225, 300]);set(gcf,'color','w');
subplot(2,1,1)
plot(Ephys(temp(cnr)).IV.traces(:,14),'Color',([27 30 187]./256),'LineWidth',1);set(gca,'box','off');axis off;
hold on;ylim([-150 50]);
%hold on;text(-200*srF,Ephys(temp(cnr)).IV.RMP,[num2str(Ephys(temp(cnr)).IV.RMP),'mV'],'FontSize',9);
subplot(2,1,2)
plot(Ephys(temp(cnr2)).IV.traces(:,10),'Color',([27 117 187]./256),'LineWidth',1);set(gca,'box','off');axis off;
hold on;ylim([-150 50]);
%hold on;text(-200*srF,Ephys(temp(cnr2)).IV.RMP,[num2str(Ephys(temp(cnr2)).IV.RMP),'mV'],'FontSize',9);


%% scatter plot (spike freq vs EPSP amplitude)
pointsize=100;
fig4=figure;set(fig4, 'Position', [200, 400, 550, 500]);set(gcf,'color','w');
% s1=scatter(cpn_spk,epsp_cpn,pointsize,'c','filled','MarkerFaceColor',([0 166 156]./256),'MarkerEdgeColor','w');hold on;
% s2=scatter(ntsr_spk,epsp_ntsr,pointsize,'m','filled','MarkerFaceColor',([235 0 139]./256),'MarkerEdgeColor','w');hold on;
% s3=scatter(penk_spk,epsp_penk,pointsize,'r','filled','MarkerFaceColor',([190 30 45]./256),'MarkerEdgeColor','w');hold on;
%s4=scatter(cci_spk,epsp_cci,pointsize,'r','filled','MarkerFaceColor',([190 30 45]./256),'MarkerEdgeColor','w');hold on;
s5=scatter(pv_spk,epsp_pv,pointsize,'y','filled','MarkerFaceColor',([246 146 30]./256),'MarkerEdgeColor','w');hold on;
s6=scatter(som_spk,epsp_som,pointsize,'b','filled','MarkerFaceColor',([27 117 187]./256),'MarkerEdgeColor','w');hold on;
s7=scatter(vip_spk,epsp_vip,pointsize,'g','filled','MarkerFaceColor',([0 165 81]./256),'MarkerEdgeColor','w');hold on;
s8=scatter(gad_spk,epsp_gad,pointsize,'k','filled','MarkerFaceColor',([0.5 0.5 0.5]),'MarkerEdgeColor','w');hold on;
xlim([-3 120]);ylim([-3 120]);
hold on;line([threshold_FS threshold_FS],[ylim],'linewidth',0.5,'Color','k','LineStyle','--');
hold on;line([xlim],[45 45],'linewidth',0.5,'Color','k','LineStyle','--');
% legend({'CPN','CT','CCi','CCi','PV','SST','VIP','GAD'},'Location','northwest');
legend({'PV','SST','VIP','GAD'},'Location','northwest');
legend boxoff
box off; ylabel('Light evoked response (mV)');xlabel('Max. spike frequency (Hz)');
set(gca,'FontSize',12);


%% Example spiking 


ov_min=-5;ov_max=20;
range=4000:8000;

%CPN
cnr=88;
temp=[];temp=find(cpn_k_cells==1);
fig4=figure;set(fig4, 'Position', [200, 200, 430, 300]);set(gcf,'color','w');
subplot(3,2,1)
plot(Ephys(temp(cnr)).sub_traces_train(range,2),'Color','k','LineWidth',1.3);set(gca,'box','off');
hold on;plot([1000 1000],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
%end
ylim([ov_min-10 ov_max]);title('CPN','Color',([0 166 156]./256));
axis off;

%CT
cnr=17;%
temp=[];temp=find(ntsr_k_cells==1);
subplot(3,2,3)
plot(Ephys(temp(cnr)).sub_traces_train(range,2),'Color','k','LineWidth',1.3);set(gca,'box','off');
%hold on;plot([1000 1000],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
%end
ylim([ov_min-10 ov_max]);title('CT','Color',([235 0 139]./256));
axis off;

%Penk
cnr=5;%
temp=[];temp=find(penk_k_cells==1);
subplot(3,2,5)
plot(Ephys(temp(cnr)).sub_traces_train(range,2),'Color','k','LineWidth',1.3);set(gca,'box','off');
%hold on;plot([1000 1000],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
ylim([ov_min-10 ov_max]);title('CCi','Color',([190 30 45]./256));
%axis off;

%PV
cnr=12;%
ov_min=-5;ov_max=100;
temp=[];temp=find(pv_k_cells==1);
subplot(3,2,2)
plot(Ephys(temp(cnr)).sub_traces_train(range,2),'Color','k','LineWidth',1.3);set(gca,'box','off');
hold on;plot([1000 1000],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
ylim([ov_min-10 ov_max]);title('PV','Color',([246 146 30]./256));
axis off;

%SOM
cnr=22;%
ov_min=-5;ov_max=100;
temp=[];temp=find(som_k_cells==1);
subplot(3,2,4)
plot(Ephys(temp(cnr)).sub_traces_train(range,2),'Color','k','LineWidth',1.3);set(gca,'box','off');
%hold on;plot([1000 1000],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
ylim([ov_min-10 ov_max]);title('SST','Color',([27 117 187]./256));
axis off;

%VIP
cnr=4;%
ov_min=-5;ov_max=100;
temp=[];temp=find(vip_k_cells==1);
subplot(3,2,6)
plot(Ephys(temp(cnr)).sub_traces_train(range,2),'Color','k','LineWidth',1.3);set(gca,'box','off');
%hold on;plot([1000 1000],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
ylim([ov_min-10 ov_max]);title('VIP','Color',([0 165 81]./256));
%axis off;
%% TTX experiemnts 

%NON PAIRED CS solution without TTX 4AP, including wash in
temp2=[];cpn_ttx=[];
temp2=cell_selecter(Ephys,'label',1,'sol',2,'drugs',1);
cpn_ttx=temp2;

temp=[];temp=find(cpn_ttx==1);ttx_ipsc=[];ttx_epsc=[];

cnr=3;
ov_min=-20;ov_max=600;
start=4000;
endp=10000;
fig4=figure;set(fig4, 'Position', [200, 400, 150, 300]);set(gcf,'color','w');
subplot(2,1,1)
plot(Ephys(temp(cnr)).sub_traces_train(start:endp,3),'Color','k','LineWidth',1.5);set(gca,'box','off');
hold on;plot(Ephys(temp(cnr)).sub_traces_train(start:endp,2),'Color','b','LineWidth',1.5);set(gca,'box','off');
hold on;plot([1000 1000],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
text(1500,-40,'TTX + 4AP');hold on;text(1500,450,'IPSC before','Color','b');
ylim([-80 600]);
%axis off; 
set(gca,'FontSize',10);
subplot(2,1,2)
 start=24000;
 endp=30000;
plot(Ephys(temp(cnr)).sub_traces_train(start:endp,4),'Color','k','LineWidth',1.5);set(gca,'box','off');
hold on;plot(Ephys(temp(cnr)).sub_traces_train(start:endp,1),'Color','b','LineWidth',1.5);set(gca,'box','off');
ylim([-150 10]);
%hold on;plot([1000 1000],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
%text(1500,-40,'TTX + 4AP');hold on;text(1500,450,'IPSC before','Color','b');
%% 

ttx_ipsc=[max(abs(Ephys(temp(1)).train_p(1:2,2))) max(abs(Ephys(temp(1)).train_p(1:2,3)));...
    max(abs(Ephys(temp(3)).train_p(1:2,2))) max(abs(Ephys(temp(3)).train_p(1:2,3)));...
  max(abs(Ephys(temp(10)).train_p(1:2,2))) max(abs(Ephys(temp(10)).train_p(1:2,3)));...
  max(abs(Ephys(temp(11)).train_p(1:2,2))) max(abs(Ephys(temp(11)).train_p(1:2,3)));...
  max(abs(Ephys(temp(12)).train_p(1:2,2))) max(abs(Ephys(temp(12)).train_p(1:2,3)));...
   max(abs(Ephys(temp(5)).train_p(1:2,1))) max(abs(Ephys(temp(5)).train_p(1:2,2)))];

% ttx_epsc=[max(abs(Ephys(temp(1)).high_n(1:2,1))) max(abs(Ephys(temp(1)).high_n(1:2,4)));...
%     max(abs(Ephys(temp(2)).high_n(1:2,1))) max(abs(Ephys(temp(2)).high_n(1:2,4)));...
%   max(abs(Ephys(temp(3)).high_n(1:2,1))) max(abs(Ephys(temp(3)).high_n(1:2,4)));...
%    max(abs(Ephys(temp(4)).high_n(1:2,1))) max(abs(Ephys(temp(4)).high_n(1:2,3)))];

cl={'b','k'};
data=[];data=ttx_ipsc;
paired_plot_box(data,cl);
xticklabels({'before','TTX + 4AP'});ylabel('IPSC amplitude (pA)');set(gca,'FontSize',10);
xtickangle(45);%yticks([0:125:250]); set(gca,'FontSize',10);
ylim([0 600]);yticks([0:200:600]);
 set(gca,'FontSize',10);
% cl={'r','k'};
% data=[];data=ttx_epsc;
% paired_plot_box(data,cl);
% xticklabels({'before','TTX + 4AP'});ylabel('EPSC amplitude (pA)');set(gca,'FontSize',10);
% xtickangle(45); set(gca,'FontSize',10);
% yticks([0:100:300]);





