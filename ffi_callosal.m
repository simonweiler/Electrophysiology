%% FEEDFORWARD INHIBITION
%which cell types spike upon callosal activation in the slice? : CPN, CT,
%CCi (Penk and unlablled CCs), PV, SST, VIP, GAD
%% Load data structure 
str_L6    = 'D:\Postdoc_Margrie\Projects\Callosal\output';
folder_list = uipickfiles('FilterSpec',str_L6);
load(char(folder_list));
%sampling rate
srF=20;
sr=20000;
%% Read out all EPSPs across cell types

stim_type=1;

%PV
pv_k_cells=cell_selecter(Ephys,'drugs',0,'label',3,'sol',1);
epsp_pv=[];
[epsp_pv] = readout_amp_epsp(Ephys,pv_k_cells,stim_type,sr);

% CPN
cpn_k_cells=cell_selecter(Ephys,'drugs',0,'label',1,'sol',1);
[epsp_cpn]=[];
[epsp_cpn] = readout_amp_epsp(Ephys,cpn_k_cells,stim_type,sr);

% Ntsr1
ntsr_k_cells=cell_selecter(Ephys,'drugs',0,'label',4,'geno',6,'sol',1);
[epsp_ntsr]=[];
[epsp_ntsr] = readout_amp_epsp(Ephys,ntsr_k_cells,stim_type,sr);

% Penk
penk_k_cells=cell_selecter(Ephys,'drugs',0,'label',5,'geno',7,'sol',1);
[epsp_penk]=[];
[epsp_penk] = readout_amp_epsp(Ephys,penk_k_cells,stim_type,sr); 

% Unlabelled excitatry cell in ntsr1 animal (additional to Penk CCi) only K solution
cci_k_cells=cell_selecter(Ephys,'drugs',0,'label',0,'geno',6,'sol',1);
[epsp_cci]=[];
[epsp_cci] = readout_amp_epsp(Ephys,cci_k_cells,stim_type,sr); 

% SOM
som_k_cells=cell_selecter(Ephys,'drugs',0,'label',6,'geno',8,'sol',1,'optovariant',1);
[epsp_som]=[];
[epsp_som] = readout_amp_epsp(Ephys,som_k_cells,stim_type,sr);

% VIP
vip_k_cells=cell_selecter(Ephys,'drugs',0,'label',8,'geno',9,'sol',1,'optovariant',1);
[epsp_vip]=[];
[epsp_vip] = readout_amp_epsp(Ephys,vip_k_cells,stim_type,sr);

% read out all GAD
in_k_cells=cell_selecter(Ephys,'drugs',0,'label',2,'sol',1);
[epsp_gad]=[];
[epsp_gad] = readout_amp_epsp(Ephys,in_k_cells,stim_type,sr);
%% Read out max spike frequency for all cell types
[cpn_spk] = spike_counts(Ephys,cpn_k_cells);
[ntsr_spk] = spike_counts(Ephys,ntsr_k_cells);
[penk_spk] = spike_counts(Ephys,penk_k_cells);
[cci_spk] = spike_counts(Ephys,cci_k_cells);
[pv_spk] = spike_counts(Ephys,pv_k_cells);
[som_spk] = spike_counts(Ephys,som_k_cells);
[vip_spk] = spike_counts(Ephys,vip_k_cells);
[gad_spk] = spike_counts(Ephys,in_k_cells);
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
hold on;line([30 30],[ylim],'linewidth',0.5,'Color','k','LineStyle','--');
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





