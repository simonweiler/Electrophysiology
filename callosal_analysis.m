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
%Ntsr1 mouse line, K-gluc, retro cells
retro_cs = cell_selecter(Ephys, 'label',1, 'sol',2,'geno',6);
%Ntsr1 mouse line, K-gluc, ntsr1+ cells
ntsr_cs = cell_selecter(Ephys, 'label',4, 'sol',2,'geno',6);

%PV mouse line, K-gluc, retro cells
retro_k_pv = cell_selecter(Ephys, 'label',1, 'sol',1,'geno',5);
%PV mouse line, K-gluc, ntsr1+ cells
pv_k = cell_selecter(Ephys, 'label',3, 'sol',1,'geno',5);
%PV mouse line, K-gluc, retro cells
retro_cs_pv = cell_selecter(Ephys, 'label',1, 'sol',2,'geno',5);
%PV mouse line, K-gluc, ntsr1+ cells
pv_cs = cell_selecter(Ephys, 'label',3, 'sol',2,'geno',5);
%% Instrinisc properties Ntsr1 vs retro
temp=[];
for i=1:length(find(retro_k==1));
    temp=find(retro_k==1);
    rmp_retro(i)=Ephys(temp(i)).IV.RMP;
end
temp=[];
for i=1:length(find(ntsr_k==1));
    temp=find(ntsr_k==1);
    rmp_ntsr(i)=Ephys(temp(i)).IV.RMP;
   
end
rmp_1= [rmp_retro' rmp_ntsr'];
%% 

fig3= figure;set(fig3, 'Name', 'RMP');set(fig3, 'Position', [200, 300, 200, 200]);set(gcf,'color','w');
hold on;
for i=1:length(rmp_1)
     p1=plot([1,2],[rmp_1(:,1),rmp_1(:,2)],'color',[0.5 0.5 0.5]);
     
end

hold on;pS=plotSpread([rmp_1(:,1),rmp_1(:,2)],'categoryIdx',[ones(1,length(rmp_1(:,1)))' ones(1,length(rmp_1(:,2)))'*2],...
    'categoryMarkers',{'o','o'},'categoryColors',{'k','r'});hold on;
hold on;plot([1,2],[nanmedian(rmp_1(:,1)),nanmedian(rmp_1(:,2))],'k','LineWidth',3);
box off;xticklabels({'CPN','NtsR1'});ylabel('RMP (mV)');yticks([-80:5:-60]);set(gca,'FontSize',10);
%set(gca, 'YScale', 'log')
[p1]=signrank(rmp_1(:,1) ,rmp_1(:,2));title([' p=' num2str(p1) ', n=' num2str(length(rmp_1))]);
