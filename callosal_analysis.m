%Callosal input mapping and ephys charaterization  scripts 
%genotypes 1=Tlx3Cre; 2=Rpb4Cre; 3=Syt6Cre; 4=GADcre tdtomato; 5= PVcre tdtomat; 6=Ntsr1cre tdtomato 
%% Load data structure 
str_L6    = 'D:\Postdoc_Margrie\Projects\Callosal\output';
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
%% nonlabelled cells
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
%% Plot IV for all reto in ntsr mouse line cells
plot_intrinsic(Ephys,retro_k,4,srF,'k',1);
%% Plot Rheobase 1x and 2x for all Ntsr1+ cells
[F2xRheo]=plot_intrinsic(Ephys,ntsr_k,4,srF,'k',2);
%% Plot Rheobase 1x and 2x for all retro+ cells
[F2xRheo_rk]=plot_intrinsic(Ephys,retro_k,4,srF,'m',2);
%% Plot Rheobase 1x and 2x for all nonlabelled cells
plot_intrinsic(Ephys,nlk,5,srF,'m',2);
%% Time fraction between first two spiked
tim_1=[F2xRheo_rk' F2xRheo'];
%% Read out intrinsic properties for all retro cells with K-gluc across all mouse lines
%extract info from IV trace and passive, extract Rheobase traces where at
%least 2 spikes are present
[rmp_rk maxsp_rk rheo_rk rin_rk tau_rk sag_rk trace_rk st_rk] = passive_readout(Ephys,rk)
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
%% Instrinisc properties Ntsr1 vs retro in the Ntsr mouse line
%retro
[rmp_retro maxsp_retro rheo_retro rin_retro tau_retro sag_retro trace_retro st_retro] = passive_readout(Ephys,retro_k);
%ntsr
[rmp_ntsr maxsp_ntsr rheo_ntsr rin_ntsr tau_ntsr sag_ntsr trace_ntsr st_ntsr] = passive_readout(Ephys,ntsr_k);
rmp_1= [rmp_retro' rmp_ntsr'];
maxsp_1=[maxsp_retro' maxsp_ntsr'];
rheo_1=[rheo_retro' rheo_ntsr'];
rin_1=[rin_retro' rin_ntsr'];
tau_1=[tau_retro' tau_ntsr'];
sag_1=[sag_retro' sag_ntsr'];
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
p1=paired_plot(rmp_1,0);xticklabels({'CPN','NtsR1'});ylabel('RMP (mV)');yticks([-80:5:-60]);set(gca,'FontSize',10);
%Input resistance Rin
p2=paired_plot(rin_1,0);xticklabels({'CPN','NtsR1'});ylabel('Input resistance (mOhm)');set(gca,'FontSize',10);
%Tau 
p3=paired_plot(tau_1,0);xticklabels({'CPN','NtsR1'});ylabel('Tau (ms)');set(gca,'FontSize',10);
%Maximum spike rate
p4=paired_plot(maxsp_1,0);xticklabels({'CPN','NtsR1'});ylabel('Max spike number');yticks([0:5:30]);set(gca,'FontSize',10);
%Rheobase
p5=paired_plot(rheo_1,1);xticklabels({'CPN','NtsR1'});ylabel('Rheobase (pA)');yticks([0:50:200]);set(gca,'FontSize',10);
%Sag Ratio
p6=paired_plot(sag_1,1);xticklabels({'CPN','NtsR1'});ylabel('Sag Ratio');set(gca,'FontSize',10);
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
%% Plot EPSC and IPSC retro vs NtsR1
%EPSC
paired_plot(epsc_1,0);xticklabels({'CPN','NtsR1'});ylabel('EPSC peak (pA)');set(gca,'FontSize',10);
%IPSC
paired_plot(ipsc_1,0);xticklabels({'CPN','NtsR1'});ylabel('IPSC peak (pA)');set(gca,'FontSize',10);
%E_I
paired_plot(e_i_1,0);xticklabels({'CPN','NtsR1'});ylabel('E / I Ratio');set(gca,'FontSize',10);
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
