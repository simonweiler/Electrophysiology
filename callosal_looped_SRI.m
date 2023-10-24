%This script reads out the peak amplitude and E/I ratios for synaptic input for different
%cell classes, 

% calculates the CRI= (EPSC non looped - EPSPC looped / EPSC non looped +
% EPSPC looped)

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
%% readout input strengths for different groups and internal solutions 
[ntsr_cs ntsr_cpncs ntsr_k ntsr_cpnk ...
    penk_cs penk_cpncs penk_k penk_cpnk ...
    cci_k cci_cpnk...
    som_cs som_cpncs som_k som_cpnk ...
    pv_cs pv_cpncs pv_k pv_cpnk vip_cs vip_cpncs vip_k vip_cpnk] = readout_input_strength_callosal(Ephys);
%% Calculate CRI ratios (CRACM reponse index)
%NTSR
ntsr_ratios=[];
for i=1:size(ntsr_cs,2)
ntsr_er=[];
ntsr_er=[(ntsr_cs(:,i)-ntsr_cpncs(:,i))./(ntsr_cs(:,i)+ntsr_cpncs(:,i)) ;...
    (ntsr_k(:,i)-ntsr_cpnk(:,i))./(ntsr_k(:,i)+ntsr_cpnk(:,i))];
ntsr_er(find(ntsr_er==Inf))=NaN;
ntsr_ratios(:,i)=ntsr_er;
end

%PENK
penk_ratios=[];
for i=1:size(penk_cs,2)
penk_er=[];
penk_er=[(penk_cs(:,i)-penk_cpncs(:,i))./(penk_cs(:,i)+penk_cpncs(:,i)) ;...
    (penk_k(:,i)-penk_cpnk(:,i))./(penk_k(:,i)+penk_cpnk(:,i))];
penk_er(find(penk_er==Inf))=NaN;
penk_ratios(:,i)=penk_er;
end

%CCI
cci_ratios=[];
for i=1:size(cci_k,2)
cci_er=[];
cci_er=[(cci_k(:,i)-cci_cpnk(length(cci_k):end,1))./(cci_k(:,i)+cci_cpnk(length(cci_k):end,1))];
cci_er(find(cci_er==Inf))=NaN;
cci_ratios(:,i)=cci_er;
end

%PV
pv_ratios=[];
for i=1:size(pv_cs,2)
pv_er=[];
pv_er=[(pv_cs(:,i)-pv_cpncs(:,i))./(pv_cs(:,i)+pv_cpncs(:,i)) ;...
    (pv_k(:,i)-pv_cpnk(:,i))./(pv_k(:,i)+pv_cpnk(:,i))];
pv_er(pv_er==Inf)=NaN;
pv_ratios(:,i)=pv_er;
end

%SOM
som_ratios=[];
for i=1:size(som_cs,2)
som_er=[];
som_er=[(som_cs(:,i)-som_cpncs(:,i))./(som_cs(:,i)+som_cpncs(:,i)) ;...
    (som_k(:,i)-som_cpnk(:,i))./(som_k(:,i)+som_cpnk(:,i))];
som_er(find(som_er==Inf))=NaN;
som_ratios(:,i)=som_er;
end
% VIP
vip_ratios=[];
for i=1:size(vip_cs,2)
vip_er=[];
vip_er=[(vip_cs(:,i)-vip_cpncs(:,i))./(vip_cs(:,i)+vip_cpncs(:,i)) ;...
    (vip_k(:,i)-vip_cpnk(:,i))./(vip_k(:,i)+vip_cpnk(:,i))];
vip_er(find(vip_er==Inf))=NaN;
vip_ratios(:,i)=vip_er;
end
%% combine PENK and CCi
cci_all=[cci_ratios ;penk_ratios];

%% check all statistically
%define which pulses
col=7;col2=8;%thrs=1000;
ct_data=[[ntsr_cpncs(:,col) ;ntsr_cpnk(:,col)] [ntsr_cs(:,col) ;ntsr_k(:,col)]];
cci_data=[[penk_cpncs(:,col) ;penk_cpnk(:,col); cci_cpnk(length(cci_k):end,col)] [penk_cs(:,col) ;penk_k(:,col); cci_k(:,col)]];
pv_data=[[pv_cpncs(:,col) ;pv_cpnk(:,col)] [pv_cs(:,col) ;pv_k(:,col)]];
som_data=[[som_cpncs(:,col) ;som_cpnk(:,col)] [som_cs(:,col) ;som_k(:,col)]];
vip_data=[[vip_cpncs(:,col) ;vip_cpnk(:,col)] [vip_cs(:,col) ;vip_k(:,col)]];

ct_data2=[[ntsr_cpncs(:,col2) ;ntsr_cpnk(:,col2)] [ntsr_cs(:,col2) ;ntsr_k(:,col2)]];
cci_data2=[[penk_cpncs(:,col2) ;penk_cpnk(:,col2); cci_cpnk(length(cci_k):end,col2)] [penk_cs(:,col2) ;penk_k(:,col2); cci_k(:,col2)]];
pv_data2=[[pv_cpncs(:,col2) ;pv_cpnk(:,col2)] [pv_cs(:,col2) ;pv_k(:,col2)]];
som_data2=[[som_cpncs(:,col2) ;som_cpnk(:,col2)] [som_cs(:,col2) ;som_k(:,col2)]];
vip_data2=[[vip_cpncs(:,col2) ;vip_cpnk(:,col2)] [vip_cs(:,col2) ;vip_k(:,col2)]];
% if sum(any(vip_data>thrs))>0
%     [m l]=find(vip_data>thrs);
%     vip_data(m,l)=vip_data2(m,l);
% end
% if sum(any(pv_data>thrs))>0
%     [m l]=find(pv_data>thrs);
%     pv_data(m,l)=pv_data2(m,l);
% end
[p_ct]=ranksum(ct_data(:,1) ,ct_data(:,2));
[p_ct]=ranksum(ct_data2(:,1) ,ct_data2(:,2));
[p_cci]=ranksum(cci_data(:,1) ,cci_data(:,2));
[p_cci2]=ranksum(cci_data2(:,1) ,cci_data2(:,2));
[p_pv]=ranksum(pv_data(:,1) ,pv_data(:,2));
[p_pv2]=ranksum(pv_data2(:,1) ,pv_data2(:,2));
[p_som]=ranksum(som_data(:,1) ,som_data(:,2));
[p_som2]=ranksum(som_data2(:,1) ,som_data2(:,2));
[p_vip]=ranksum(vip_data(:,1) ,vip_data(:,2));
[p_vip2]=ranksum(vip_data2(:,1) ,vip_data2(:,2));
%% New aprroach: find the the maximum for 1 Hz, 10 Hz or 25 Hz and for looped CPN cell and then find corresponding pulse for respective non-looped cells 
%col: 1 4 7 are excitatory 1, 10 or 25 Hz
%col: 2 5 8 are inhbitory 1, 10 or 25 Hz

%col=[1 4 7];
col=[2 5 8];

%Ntsr1 aka CT
[ct_m ct_i]=max([ntsr_cpncs(:,col) ;ntsr_cpnk(:,col)],[],2);
%[ct_m ct_i]=max([ntsr_cs(:,col) ;ntsr_k(:,col)],[],2);
ct_temp=[ntsr_cs(:,col) ;ntsr_k(:,col)];
%ct_temp=[ntsr_cpncs(:,col) ;ntsr_cpnk(:,col)];
ct_r=[ntsr_ratios(:,col)];
ct_m2=[];
ct_cri=[];
for i=1:length(ct_temp)
    %amplitude
    ct_m2(i)=ct_temp(i,ct_i(i));
     %ratio
    ct_cri(i)=ct_r(i,ct_i(i));
end
ct_data=[];
ct_data=[ct_m ct_m2'];
%% 

%Penk and non CT non CPN
%[cc_m cc_i]=max([penk_cpncs(:,col) ;penk_cpnk(:,col); cci_cpnk(length(cci_k):end,col)],[],2);
[cc_m cc_i]=max([penk_cpncs(:,col) ;penk_cpnk(:,col); cci_cpnk(:,col)],[],2);
%[cc_m cc_i]=max([penk_cs(:,col) ;penk_k(:,col); cci_k(length(cci_k):end,col)],[],2);
cc_temp=[penk_cs(:,col) ;penk_k(:,col); cci_k(:,col)];
%cc_temp=[penk_cpncs(:,col) ;penk_cpnk(:,col); cci_cpnk(:,col)];
cc_r=[cci_all(:,col)];
cc_m2=[];
cc_cri=[];
for i=1:length(cc_i)
    cc_m2(i)=cc_temp(i,cc_i(i));
    cc_cri(i)=cc_r(i,cc_i(i));
end
cci_data=[];
cci_data=[cc_m cc_m2'];
%% 

%PV
[pv_m pv_i]=max([pv_cpncs(:,col) ;pv_cpnk(:,col)],[],2);
pv_temp=[pv_cs(:,col) ;pv_k(:,col)];
pv_r=[pv_ratios(:,col)];
pv_m2=[];
pv_cri=[];
for i=1:length(pv_temp)
    pv_m2(i)=pv_temp(i,pv_i(i));
    pv_cri(i)=pv_r(i,pv_i(i));
end
pv_data=[];
pv_data=[pv_m pv_m2'];
%% 

%SOM
[som_m som_i]=max([som_cpncs(:,col) ;som_cpnk(:,col)],[],2);
som_temp=[som_cs(:,col) ;som_k(:,col)];
som_r=[som_ratios(:,col)];
som_m2=[];
som_cri=[];
for i=1:length(som_temp)
    som_m2(i)=som_temp(i,som_i(i));
    som_cri(i)=som_r(i,som_i(i));
end
som_data=[];
som_data=[som_m som_m2'];
%% 

%VIP
[vip_m vip_i]=max([vip_cpncs(:,col) ;vip_cpnk(:,col)],[],2);
vip_temp=[vip_cs(:,col) ;vip_k(:,col)];
vip_r=[vip_ratios(:,col)];
vip_m2=[];
vip_cri=[];
for i=1:length(vip_temp)
    %amplitude
    vip_m2(i)=vip_temp(i,vip_i(i));
    %ratio
    vip_cri(i)=vip_r(i,vip_i(i));
end
vip_data=[];
vip_data=[vip_m vip_m2'];
%% Do statistics ranksum 
[p_ct]=ranksum(ct_data(:,1) ,ct_data(:,2))
[p_cci]=ranksum(cci_data(:,1) ,cci_data(:,2))
[p_pv]=ranksum(pv_data(:,1) ,pv_data(:,2))
[p_som]=ranksum(som_data(:,1) ,som_data(:,2))
[p_vip]=ranksum(vip_data(:,1) ,vip_data(:,2))
%% OTHER OPTION USE MEAN ACROSS 1, 5 and 10 Hz
col=[2 5 8];
ct_avg3=[];cc_avg3=[];pv_avg3=[];som_avg3=[];vip_avg3=[];
ct_avg3=[nanmean([ntsr_cpncs(:,col) ;ntsr_cpnk(:,col)],2) nanmean([ntsr_cs(:,col) ;ntsr_k(:,col)],2)];
[p_ct]=ranksum(ct_avg3(:,1) ,ct_avg3(:,2))
cc_avg3=[nanmean([penk_cpncs(:,col) ;penk_cpnk(:,col); cci_cpnk(:,col)],2) nanmean([penk_cs(:,col) ;penk_k(:,col); cci_k(:,col)],2)]
[p_cc]=ranksum(cc_avg3(:,1) ,cc_avg3(:,2))
pv_avg3=[nanmean([pv_cpncs(:,col) ;pv_cpnk(:,col)],2) nanmean([pv_cs(:,col) ;pv_k(:,col)],2)];
[p_pv]=ranksum(pv_avg3(:,1) ,pv_avg3(:,2))
som_avg3=[nanmean([som_cpncs(:,col) ;som_cpnk(:,col)],2) nanmean([som_cs(:,col) ;som_k(:,col)],2)];
[p_som]=ranksum(som_avg3(:,1) ,som_avg3(:,2))
vip_avg3=[nanmean([vip_cpncs(:,col) ;vip_cpnk(:,col)],2) nanmean([vip_cs(:,col) ;vip_k(:,col)],2)];
[p_vip]=ranksum(vip_avg3(:,1) ,vip_avg3(:,2))
%% 
ct_cri=[];ct_cri=(ct_avg3(:,2)-ct_avg3(:,1))./(ct_avg3(:,2)+ct_avg3(:,1));
cc_cri=[];cc_cri=(cc_avg3(:,2)-cc_avg3(:,1))./(cc_avg3(:,2)+cc_avg3(:,1));
pv_cri=[];pv_cri=(pv_avg3(:,2)-pv_avg3(:,1))./(pv_avg3(:,2)+pv_avg3(:,1));
som_cri=[];som_cri=(som_avg3(:,2)-som_avg3(:,1))./(som_avg3(:,2)+som_avg3(:,1));
vip_cri=[];vip_cri=(vip_avg3(:,2)-vip_avg3(:,1))./(vip_avg3(:,2)+vip_avg3(:,1));
%% Violin plots for CRI
g1=[];g2=[];g3=[];g4=[];g5=[];
p1=[];p2=[];p3=[];p4=[];p5=[];
p1=ct_cri;p2=cc_cri;p3=pv_cri;p4=som_cri;p5=vip_cri;
par=[];par=[p1 p2 p3 p4 p5]';
g1=[];g1=ones(1,length(p1));
g2=[];g2=ones(1,length(p2))*2;
g3=[];g3=ones(1,length(p3))*3;
g4=[];g4=ones(1,length(p4))*4;
g5=[];g5=ones(1,length(p5))*5;
gro=[];gro=[g1 g2 g3 g4 g5]';
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 300, 300]);set(gcf,'color','w');
hold on;line([0 6],[0 0],'LineStyle',':','Color','k','LineWidth',1)
violins = violinplot(par, gro,'ViolinColor',[([235 0 139]./256);([190 30 45]./256);([246 146 30]./256);([27 117 187]./256);([0 165 81]./256)]);box off;
xlim([0 6]);ylabel('CRI');set(gca,'FontSize',12);
% h = gca;h.XAxis.Visible = 'off';
% set(gca,'xtick',[]);
xticklabels({'CT','CCi','PV','SST','VIP'});
ylim([-1.2 1.2]);
%% CT as example in paried scatter plot
fig1= figure;set(fig1, 'Name', 'Barplot groups');set(fig1, 'Position', [400, 200, 300, 300]);set(gcf,'color','w');
% scatter(ct_data(:,1)./max(ct_data(:,1)),ct_data(:,2)./max(ct_data(:,2)),20,'b','filled','o');axis square;hold on;rf=refline(1,0); rf.Color= [0 0 0];rf.LineStyle= '--';
% hold on;scatter(cci_data(:,1)./max(cci_data(:,1)),cci_data(:,2)./max(cci_data(:,2)),20,'r','filled','o');axis square;hold on;rf=refline(1,0); rf.Color= [0 0 0];rf.LineStyle= '--';
% hold on;scatter(pv_data(:,1)./max(pv_data(:,1)),pv_data(:,2)./max(pv_data(:,2)),20,'y','filled','o');axis square;hold on;rf=refline(1,0); rf.Color= [0 0 0];rf.LineStyle= '--';
% hold on;scatter(som_data(:,1)./max(som_data(:,1)),som_data(:,2)./max(som_data(:,2)),20,'m','filled','o');axis square;hold on;rf=refline(1,0); rf.Color= [0 0 0];rf.LineStyle= '--';
% hold on;scatter(vip_data(:,1)./max(vip_data(:,1)),vip_data(:,2)./max(vip_data(:,2)),20,'g','filled','o');axis square;hold on;rf=refline(1,0); rf.Color= [0 0 0];rf.LineStyle= '--';
 s1=scatter(ct_data(:,1),ct_data(:,2),35,'b','filled','o');axis square;hold on;rf=refline(1,0); rf.Color= [0 0 0];rf.LineStyle= '--';
 s1.MarkerEdgeColor = [0 0 0];
 %xlim([0 400]);ylim([0 400]);
 xlabel('Light-evoked amplitude CPN (pA)');ylabel('Light-evoked amplitude CT (pA)');
text(50,900,['p=' num2str(round(p_ct,2))]);
hold on;plot([0 max([xlim ylim])], [0 max([xlim ylim])], '--k');
set(gca,'FontSize',12);
 % hold on;scatter(cci_data(:,1),cci_data(:,2),20,'g','filled','o');axis square;hold on;rf=refline(1,0); rf.Color= [0 0 0];rf.LineStyle= '--';
%   hold on;scatter(pv_data(:,1),pv_data(:,2),20,'b','filled','o');axis square;hold on;rf=refline(1,0); rf.Color= [0 0 0];rf.LineStyle= '--';
% hold on;scatter(som_data(:,1),som_data(:,2),20,'y','filled','o');axis square;hold on;rf=refline(1,0); rf.Color= [0 0 0];rf.LineStyle= '--';
%    hold on;scatter(vip_data(:,1),vip_data(:,2),20,'g','filled','o');axis square;hold on;rf=refline(1,0); rf.Color= [0 0 0];rf.LineStyle= '--';
%    hold on;scatter(log(pv_data(:,1)),log(pv_data(:,2)),20,'b','filled','o');axis square;hold on;rf=refline(1,0); rf.Color= [0 0 0];rf.LineStyle= '--';
%   hold on;scatter(log(som_data(:,1)),log(som_data(:,2)),20,'y','filled','o');axis square;hold on;rf=refline(1,0); rf.Color= [0 0 0];rf.LineStyle= '--';
%    hold on;scatter(log(vip_data(:,1)),log(vip_data(:,2)),20,'g','filled','o');axis square;hold on;rf=refline(1,0); rf.Color= [0 0 0];rf.LineStyle= '--';
% %% Scatter plots CRI for all 
% typ_psc=8;
% yl='E / I ratio';
% cl={'b','r','y',[1 0.5 0],'g'};
% con_data=[];con_categ=[];
% data={};data={ntsr_ratios(:,typ_psc) cci_all(:,typ_psc) pv_ratios(:,typ_psc) som_ratios(:,typ_psc) vip_ratios(:,typ_psc)};
% for i=1:length(data)
% categ{:,i}=[ones(1,length(data{:,i}))*i]';
% end
% con_data=vertcat(data{:});
% con_categ=vertcat(categ{:});
% %figure
% fig1= figure;set(fig1, 'Name', 'Barplot groups');set(fig1, 'Position', [400, 500, 290, 250]);set(gcf,'color','w');
% hold on;
%  for i=1:length(data)
%  hold on;
%  b2=plot(i,nanmean(data{:,i}),'+','MarkerSize',15);
%  hold on;
%  end
% hold on;
% plotSpread(data,'categoryIdx',con_categ,'categoryMarkers',{'o','o','o','o','o'},...
% 'categoryColors',{[0.7,0.7,0.7],[0.7,0.7,0.7],[0.7,0.7,0.7],[0.7,0.7,0.7],[0.7,0.7,0.7]});
% hold on;
% for i=1:length(data)
% er=errorbar(i,nanmean(data{:,i}),nanstd(data{:,i})/sqrt(size(data{:,i},1)));
% er.Color = [0 0 0];er.LineWidth=1;er.LineStyle = 'none'; hold on;
% end
% 
% xlim([-0 6]);
% ylim([-1.2 1.2]);
% plot(0:i+1,zeros(1,size(data,2)+2),'LineStyle',':','Color','k','LineWidth',1);
% xticks(1:5);
% xticklabels({'CT','CCi','PV','SST','VIP'});
% ylabel('CRI');
% set(gca,'FontSize',10);
%% E I ratio for E-I Index
col=1
ei_ct=(ntsr_cs(:,col)-ntsr_cs(:,col+1))./(ntsr_cs(:,col)+ntsr_cs(:,col+1));
ei_cci=(penk_cs(:,col)-penk_cs(:,col+1))./(penk_cs(:,col)+penk_cs(:,col+1));
ei_pv=(pv_cs(:,col)-pv_cs(:,col+1))./(pv_cs(:,col)+pv_cs(:,col+1));
ei_som=(som_cs(:,col)-som_cs(:,col+1))./(som_cs(:,col)+som_cs(:,col+1));
ei_vip=(vip_cs(:,col)-vip_cs(:,col+1))./(vip_cs(:,col)+vip_cs(:,col+1));
ei_cp=[(ntsr_cpncs(:,col)-ntsr_cpncs(:,col+1))./(ntsr_cpncs(:,col)+ntsr_cpncs(:,col+1));...
    (penk_cpncs(:,col)-penk_cpncs(:,col+1))./(penk_cpncs(:,col)+penk_cpncs(:,col+1));...
    (pv_cpncs(:,col)-pv_cpncs(:,col+1))./(pv_cpncs(:,col)+pv_cpncs(:,col+1));...
    (som_cpncs(:,col)-som_cpncs(:,col+1))./(som_cpncs(:,col)+som_cpncs(:,col+1))];
%% Violin plot for E-I index
g1=[];g2=[];g3=[];g4=[];g5=[];g6=[];
p1=[];p2=[];p3=[];p4=[];p5=[];p6=[];
p1=ei_cp;p2=ei_ct;p3=ei_cci;p4=ei_pv;p5=ei_som;p6=ei_vip;
par=[];par=[p1' p2' p3' p4' p5' p6']';
g1=[];g1=ones(1,length(p1));
g2=[];g2=ones(1,length(p2))*2;
g3=[];g3=ones(1,length(p3))*3;
g4=[];g4=ones(1,length(p4))*4;
g5=[];g5=ones(1,length(p5))*5;
g6=[];g6=ones(1,length(p5))*6;
gro=[];gro=[g1 g2 g3 g4 g5 g6]';
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 300, 300]);set(gcf,'color','w');
hold on;line([0 7],[0 0],'LineStyle',':','Color','k','LineWidth',1)
violins = violinplot(par, gro,'ViolinColor',[([0 166 156]./256);([235 0 139]./256);([190 30 45]./256);([246 146 30]./256);([27 117 187]./256);([0 165 81]./256)]);box off;
xlim([0 7]);ylabel('EX-IN index');set(gca,'FontSize',12);
% h = gca;h.XAxis.Visible = 'off';
% set(gca,'xtick',[]);
xticklabels({'CCc','CT','CCi','PV','SST','VIP'});
ylim([-1.2 1.2]);
%% Violin plot for excitatory cells types vs inhbitory ones
g1=[];g2=[];
p1=[];p2=[];
%combine ex and in
p1=[ei_cp ;ei_ct ;ei_cci];p2=[ei_pv ;ei_som ;ei_vip];
par=[];par=[p1' p2']';
g1=[];g1=ones(1,length(p1)); 
g2=[];g2=ones(1,length(p2))*2;
gro=[];gro=[g1 g2]';
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 300, 300]);set(gcf,'color','w');
hold on;line([0 3],[0 0],'LineStyle',':','Color','k','LineWidth',1)
violins = violinplot(par, gro,'ViolinColor',[[0.3 0.3 0.3]; [0.6 0.6 0.6]]);box off;
xlim([0 3]);ylabel('EX-IN index');set(gca,'FontSize',12);
% h = gca;h.XAxis.Visible = 'off';
% set(gca,'xtick',[]);
xticklabels({'L6 PNs','L6 INs'});
ylim([-1.2 1.2]);
%statitics
p_ei=ranksum([ei_cp ;ei_ct ;ei_cci],[ei_pv ;ei_som ;ei_vip]);











%% 
typ_psc=1;
yl='E / I ratio';
cl={'c',[0 0 1],[1 0 0],'y','m'};

data={};data={ei_cp(:,typ_psc) ei_ct(:,typ_psc) ei_cci(:,typ_psc) ei_pv(:,typ_psc) ei_som(:,typ_psc)};
for i=1:length(data)
categ{:,i}=[ones(1,length(data{:,i}))*i]';
end
con_data=vertcat(data{:});
con_categ=vertcat(categ{:});
%figure
fig1= figure;set(fig1, 'Name', 'Barplot groups');set(fig1, 'Position', [400, 500, 290, 250]);set(gcf,'color','w');
hold on;
 for i=1:length(data)
 hold on;
 b2=plot(i,nanmean(data{:,i}),'+','MarkerSize',15,'Color',cl{i});
 hold on;
 end
hold on;
plotSpread(data,'categoryIdx',con_categ,'categoryMarkers',{'o','o','o','o','o'},...
'categoryColors',{[0.7,0.7,0.7],[0.7,0.7,0.7],[0.7,0.7,0.7],[0.7,0.7,0.7],[0.7,0.7,0.7]});
hold on;
for i=1:length(data)
er=errorbar(i,nanmean(data{:,i}),nanstd(data{:,i})/sqrt(size(data{:,i},1)));
er.Color = [0 0 0];er.LineWidth=1;er.LineStyle = 'none'; hold on;
end

xlim([-0 6]);
ylim([-1.2 1.2]);
plot(0:i+1,zeros(1,size(data,2)+2),'LineStyle',':','Color','k','LineWidth',1);
xticks(1:5);
xticklabels({'CCc','CT','CCi','PV','SST'});
ylabel('EX-IN index');
set(gca,'FontSize',10);
%% AUDp VISp
ntsr_ac=cell_selecter(Ephys,'label',4,'geno',6,'sol',2,'optovariant',1,'area',3);
cp_ac=cell_selecter(Ephys,'label',1,'geno',6,'sol',2,'optovariant',1,'area',3);
%% 
[ntsr_psc] = readout_amp_update(Ephys,ntsr_ac ,1,2,1,2);
[ntsr_psc_cshf] = readout_amp_update(Ephys,ntsr_ac ,2,2,1,2);
[ntsr_psc_cshf2] = readout_amp_update(Ephys,ntsr_ac ,3,2,1,2);
[cp_psc] = readout_amp_update(Ephys,cp_ac ,1,2,1,2);
[cp_psc_cshf] = readout_amp_update(Ephys,cp_ac ,2,2,1,2);
[cp_psc_cshf2] = readout_amp_update(Ephys,cp_ac ,3,2,1,2);

%% 
col=1;
audct_data=[[cp_psc(:,col)] [ntsr_psc(:,col)]];

fig1= figure;set(fig1, 'Name', 'Barplot groups');set(fig1, 'Position', [400, 200, 300, 300]);set(gcf,'color','w');
scatter(audct_data(:,1),audct_data(:,2),20,'k','filled','o');axis square;hold on;rf=refline(1,0); rf.Color= [0 0 0];rf.LineStyle= '--';
%xlim([0 450]);ylim([0 450]);
xlabel('Looped AUDp CCc (pA)');ylabel('Non-Looped AUDp CT (pA)');


%% Scatter plots CRI for all 

%CT AUDp
audpct_ratios=[];
for i=1:size(ntsr_psc,2)
audpct_er=[];
audpct_er=[(cp_psc(:,i)-ntsr_psc(:,1))./(cp_psc(:,i)+ntsr_psc(:,1))];
audpct_er(find(audpct_er==Inf))=NaN;
audpct_ratios(:,i)=audpct_er;
end

typ_psc=1;
yl='E / I ratio';
cl={'b','r','y',[1 0.5 0],'k'};

data={};data={ntsr_ratios(:,typ_psc) cci_all(:,typ_psc) pv_ratios(:,typ_psc) som_ratios(:,typ_psc) audpct_ratios(:,typ_psc)};
for i=1:length(data)
categ{:,i}=[ones(1,length(data{:,i}))*i]';
end
con_data=vertcat(data{:});
con_categ=vertcat(categ{:});
%figure
fig1= figure;set(fig1, 'Name', 'Barplot groups');set(fig1, 'Position', [400, 500, 290, 250]);set(gcf,'color','w');
hold on;
 for i=1:length(data)
 hold on;
 b2=plot(i,nanmean(data{:,i}),'+','MarkerSize',15);
 hold on;
 end
hold on;
plotSpread(data,'categoryIdx',con_categ,'categoryMarkers',{'o','o','o','o','o'},...
'categoryColors',{[0.7,0.7,0.7],[0.7,0.7,0.7],[0.7,0.7,0.7],[0.7,0.7,0.7],[0.7,0.7,0.7]});
hold on;
for i=1:length(data)
er=errorbar(i,nanmean(data{:,i}),nanstd(data{:,i})/sqrt(size(data{:,i},1)));
er.Color = [0 0 0];er.LineWidth=1;er.LineStyle = 'none'; hold on;
end

xlim([-0 6]);
ylim([-1 1]);
plot(0:i+1,zeros(1,size(data,2)+2),'LineStyle',':','Color','k','LineWidth',1);
xticks(1:5);
xticklabels({'VISp CT','VISp CCi','VISp PV','VISp SST','AUDp CT'});
ylabel('CRI');
set(gca,'FontSize',10);