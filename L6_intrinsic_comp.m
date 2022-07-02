%callosal intrinsic electrical  properties L6 callosal cells
%% Load data structure 
str_L6    = 'D:\Postdoc_Margrie\Projects\Callosal\output';
folder_list = uipickfiles('FilterSpec',str_L6);
load(char(folder_list));
%sampling rate
srF=20;
sr=20000;
%% read out all Callosal cells with intrinsic properties across animal genotypes
cpn_k = cell_selecter(Ephys, 'label',1,'sol',1);
%% Passive properties 
[rmp_cpn maxsp_cpn rheo_cpn rin_cpn tau_cpn sag_cpn trace_cpn spike_time_cpn] = passive_readout(Ephys,cpn_k);

%% 
%remove super hyperpolarizedd cells
idx_rmv=find(rmp_cpn<-90);
rmp_cpn(idx_rmv)=NaN;maxsp_cpn(idx_rmv)=NaN;rheo_cpn(idx_rmv)=NaN;rin_cpn(idx_rmv)=NaN;tau_cpn(idx_rmv)=NaN;
sag_cpn(idx_rmv)=NaN;trace_cpn(:,idx_rmv)=NaN;spike_time_cpn(idx_rmv)=NaN;
data_pas=[];data_pas=[rmp_cpn; tau_cpn; rin_cpn; rheo_cpn;]';
stri={'V_{rest} (mV)','Tau (ms)','R_{in} (MO)','Rheobase'}
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 600, 200]);
 for i=1:size(data_pas,2)
hold on;
subplot(1,4,i)
h=histogram(data_pas(:,i),10,'FaceColor',[0.6 0.6 0.6],'EdgeColor','w','LineWidth',1)
h.EdgeColor = 'w';
h.FaceColor = [0.6 0.6 0.6];
xlabel(stri(i));
%ylim([0 1]);
%xlim([0 max(data_pass(:,i))+max(data_pass(:,i))*0.25]);
hAxis = gca;
hAxis.YAxisLocation = 'left';    % 'left' (default) or 'right'
hAxis.XAxisLocation = 'bottom'
box off
 end
 %% Active properties
spike_nr=1;
active_cpn=sp_parameters_pandora(trace_cpn,spike_nr);
%% Read out max spike frequency
 temp=[];maxsF_cpn=[];
temp=find(cpn_k==1);
spike_counts_all=[];
for i=1:length(temp)
    try
    maxsF_cpn(i)=max(Ephys(temp(i)).IV.spikecount);
    spikecount_cpn{:,i}=Ephys(temp(i)).IV.spikecount;
    temp3=Ephys(temp(i)).IV.spikecount;
    spike_countsm=ones(1,25)*NaN;
    if length(temp3)==15
    spike_countsm(1:15)=temp3;
    else
    spike_countsm(1:25)=temp3;
    end
    stimvec_cpn{:,i}=Ephys(temp(i)).IV.stimvec;
    spike_counts_all(:,i)=spike_countsm;

    catch
    maxsF_cpn(i)=NaN;
    spikecount_cpn{:,i}=NaN;
    stimvec_cpn{:,i}=NaN;
    spike_counts_all(:,i)=ones(1,25)*NaN;
    end
end
%% 
figure;imagesc((spike_counts_all./max(spike_counts_all))');xlim([7 15]);
%% 

 fig4=figure;set(fig4, 'Position', [200, 600, 300, 300]);set(gcf,'color','w');
 for i=1:length(spikecount_cpn)
     p1=plot(stimvec_cpn{:,i},spikecount_cpn{:,i},'-o');p1.Color=[0.8500 0.3250 0.0980];p1.Color(4)=3/8; 
     hold on;set(gca,'box','off');
 end
ylabel('Spike frequency (Hz)');xlabel('Injected current (pA)')
xlim([0 500]);
%% 
data_act=[];
data_act=[active_cpn([1 2 3 5 6 7 8 15],:)' maxsF_cpn'];
data_act(:,4)=data_act(:,4)*2;
data_act(:,8)=data_act(:,8)/2;
stri={'APV_{min}(mV)','APV_{peak}(mV)','APV_{thresh}(mV)', 'APVslope_{max} (\DeltamV/\Deltams)','APV_{half} (mV)','APV_{amplitude} (mV)',...
    'AHP_{max}(mV)','AP_{half width} (ms)' 'Max Spike Frequency (Hz)'};
fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 600, 400]);
 for i=1:size(data_act,2)
hold on;
subplot(3,3,i)
h=histogram(data_act(:,i),10,'FaceColor',[0.6 0.6 0.6],'EdgeColor','w','LineWidth',1)
h.EdgeColor = 'w';
h.FaceColor = [0.6 0.6 0.6];
xlabel(stri(i));
%ylim([0 1]);
%xlim([0 max(data_pass(:,i))+max(data_pass(:,i))*0.25]);
hAxis = gca;
hAxis.YAxisLocation = 'left';    % 'left' (default) or 'right'
hAxis.XAxisLocation = 'bottom'
box off
 end
%% TEST differences using PCA only including passive properties 
data_tmp=[];
data_tmp=data_pas;
data_tmp(find(sum(isnan(data_tmp),2)>0),:)=[];
score_pas=[];
[coeff_pas,score_pas,latent_pas,~,explained_pas,mu] = pca(zscore(data_tmp));
var_exp(explained_pas,[],[]); 
figure;scatter(score_pas(:,1),score_pas(:,2));
%% TEST differences using PCA only including active properties
data_tmp=[];
data_tmp=data_act;
data_tmp(find(sum(isnan(data_tmp),2)>0),:)=[];
score_act=[];
[coeff_act,score_act,latent_act,~,explained_act,mu] = pca(zscore(data_tmp));
var_exp(explained_act,[],[]); 
%figure;scatter(score_act(:,1),score_act(:,2));
score_tmp=[];score_tmp=score_act;
 fig1=figure;set(gcf,'color','w');set(fig1, 'Position', [200, 200, 250, 250])
 scatter(score_tmp(:,1),score_tmp(:,2),20,'filled', 'MarkerFaceAlpha',3/8,'MarkerFaceColor','k');
 ylabel('PC2supra');xlabel('PC1supra');set(gca,'FontSize',10);
 dip_test_SW(score_tmp(:,[1 2 3]),0,{'PC1','PC2','PC3'});ylim([0 30]);
%% Plot three examples for CPNs
 temp=[];
 step=[];step=3;
ex1=10;ex2=20;ex3=30;
temp=find(cpn_k==1);
fig4=figure;set(fig4, 'Position', [200, 800, 600, 200]);set(gcf,'color','w');
subplot(1,3,1)
plot(Ephys(temp(ex1)).IV.traces(:,1:step:end),'Color','k','LineWidth',1);set(gca,'box','off');axis off;
hold on;ylim([-150 50]);
hold on;text(-550*srF,Ephys(temp(ex1)).IV.RMP,[num2str(Ephys(temp(ex1)).IV.RMP),'mV'],'FontSize',9);
subplot(1,3,2)
plot(Ephys(temp(ex2)).IV.traces(:,1:step:end),'Color','k','LineWidth',1);set(gca,'box','off');axis off;
hold on;ylim([-150 50]);
hold on;text(-550*srF,Ephys(temp(ex2)).IV.RMP,[num2str(Ephys(temp(ex2)).IV.RMP),'mV'],'FontSize',9);
subplot(1,3,3)
plot(Ephys(temp(ex3)).IV.traces(:,1:step:end),'Color','k','LineWidth',1);set(gca,'box','off');axis off;
hold on;ylim([-150 50]);
hold on;text(-550*srF,Ephys(temp(ex3)).IV.RMP,[num2str(Ephys(temp(ex3)).IV.RMP),'mV'],'FontSize',9);
%% 
%% read out NTSR1 Cell
ntsr_k = cell_selecter(Ephys, 'label',4,'geno',6,'sol',1);
cpn_ntsr_k = cell_selecter(Ephys, 'label',1,'geno',6,'sol',1);
[rmp_ntsr maxsp_ntsr rheo_ntsr rin_ntsr tau_ntsr sag_ntsr trace_ntsr spike_time_ntsr] = passive_readout(Ephys,ntsr_k);
[rmp_ntsr_cpn maxsp_ntsr_cpn rheo_ntsr_cpn rin_ntsr_cpn tau_ntsr_cpn sag_ntsr_cpn trace_ntsr_cpn spike_time_ntsr_cpn] = passive_readout(Ephys,cpn_ntsr_k);
%% comparison 
%Input resistance 
p1=sag_ntsr_cpn;p2=sag_ntsr;
par=[];par=[p1 p2]';
g1=1:length(p1);
g2=length(p1)+1:length(par);
[statsout_rin]=dual_boxplot(par,g1,g2,0);
ylabel('Input resistance (mOhm)');set(gca,'FontSize',10);
xlim([0 3]);xticks([1 2]);xticklabels({'CPN','Ntsr1'});xtickangle(45);
%% 
spike_nr=1;
active_ntsr_cpn=sp_parameters_pandora(trace_ntsr_cpn,spike_nr);
spike_nr=1;
active_ntsr=sp_parameters_pandora(trace_ntsr,spike_nr);

%% Subplots for active paramaters 
par1=active_ntsr_cpn([1 2 3 5 6 7 8 15],:);par2=active_ntsr([1 2 3 5 6 7 8 15],:);
color_id={[0 0 0],[0.8500 0.3250 0.0980]};
str={'APV_{min}(mV)','APV_{peak}(mV)','APV_{thresh}(mV)', 'APVslope_{max} (\DeltamV/\Deltams)','APV_{half} (mV)','APV_{amplitude} (mV)',...
    'AHP_{max}(mV)','AP_{half width} (ms)'};

fig3= figure;set(fig3, 'Name', 'Paired comp');set(fig3, 'Position', [200, 500, 400, 500]);set(gcf,'color','w');
for i=1:size(par1,1)
data={};
data={par1(i,:)' par2(i,:)'};
% gr_m1=nanmean(data{:,1}); 
% gr_m2=nanmean(data{:,2}); 
subplot(3,3,i);hold on;
for k=1:length(data)
    p1=[];p2=[];
    p1=scatter(ones(1,length(data{:,1}))',data{:,1},10,'MarkerEdgeColor',color_id{1});    
    p2=scatter(ones(1,length(data{:,2}))'*2,data{:,2},10,'MarkerEdgeColor',color_id{2}); 
end
hold on;
for m=1:2
er=errorbar(m,nanmean(data{:,m}),nanstd(data{:,m})/sqrt(sum(~isnan(data{:,m}))));er.LineWidth=1.5;er.LineStyle = 'none'; hold on;
er.Color = color_id{m};er.LineWidth=1.5;er.LineStyle = 'none'; 
er.CapSize = 10;hold on;
hold on;plot(m,nanmean(data{:,m}),'o','MarkerFaceColor',color_id{m},'MarkerEdgeColor',color_id{m},'MarkerSize',7);
end
xlim([0 3]);xticks([1 2]);
xticklabels({'PN','FS'});xtickangle(45);
ylabel(str{i});
set(gca,'FontSize',10)
[p k]=ranksum(data{:,1},data{:,2});
    statsout(i)=p;
end 
%% Read out max spike frequency

