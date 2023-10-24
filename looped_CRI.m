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
srF=20;sr=20000;
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

% combine PENK and CCi
cci_all=[cci_ratios ;penk_ratios];
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
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 350, 300]);set(gcf,'color','w');
hold on;line([0 6],[0 0],'LineStyle',':','Color','k','LineWidth',1)
violins = violinplot(par, gro,'ViolinColor',[([235 0 139]./256);([190 30 45]./256);([246 146 30]./256);([27 117 187]./256);([0 165 81]./256)],'ShowMean', true,'ShowMedian', true);box off;
xlim([0 6]);ylabel('CRI');set(gca,'FontSize',12);
% h = gca;h.XAxis.Visible = 'off';
% set(gca,'xtick',[]);
xticklabels({'CT','CCi','PV','SST','VIP'});
ylim([-1.2 1.2]);
set(gca,'TickDir','out'); 
%% ex stats
%ex stats
hold on;text(4.75,1,['***'],'FontSize',18);
title('Excitation','Color','r','FontWeight','normal');
%% in stats

hold on;text(4.9,1,['*'],'FontSize',18);
hold on;text(3.9,1,['*'],'FontSize',18);
title('Inhibition','Color','b','FontWeight','normal');
%% 
%% CT as example in pairoed scatter plot excitatry
fig1= figure;set(fig1, 'Name', 'Barplot groups');set(fig1, 'Position', [400, 200, 320, 320]);set(gcf,'color','w');
 s1=scatter(ct_data(:,1),ct_data(:,2),35,'r','filled','o');axis square;hold on;rf=refline(1,0); rf.Color= [0 0 0];rf.LineStyle= '--';
 s1.MarkerEdgeColor = [0 0 0];axis square;
 xlim([0 500]);ylim([0 500]);xticks([0:250:500]);yticks([0:250:500])
 xlabel('Light-evoked amplitude CPN (pA)');ylabel('Light-evoked amplitude CT (pA)');
text(50,500,['p=' num2str(round(p_ct,2))]);
hold on;plot([0 max([xlim ylim])], [0 max([xlim ylim])], '--k');
set(gca,'FontSize',11);set(gca,'TickDir','out'); 
%% CT as example in paried scatter plot inhbitory
fig1= figure;set(fig1, 'Name', 'Barplot groups');set(fig1, 'Position', [400, 200, 300, 300]);set(gcf,'color','w');
 s1=scatter(ct_data(:,1),ct_data(:,2),35,'b','filled','o');axis square;hold on;rf=refline(1,0); rf.Color= [0 0 0];rf.LineStyle= '--';
 s1.MarkerEdgeColor = [0 0 0];
 xlim([0 1100]);ylim([0 1100]);xticks([0:500:1000]);yticks([0:500:1000])
 xlabel('Light-evoked amplitude CPN (pA)');ylabel('Light-evoked amplitude CT (pA)');
text(50,1000,['p=' num2str(round(p_ct,2))]);
hold on;plot([0 max([xlim ylim])], [0 max([xlim ylim])], '--k');
set(gca,'FontSize',11);set(gca,'TickDir','out'); title('Inhibition','Color','b','FontWeight','normal');
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
%% 
%% Violin plot for excitatory cells types vs inhbitory ones
g1=[];g2=[];
p1=[];p2=[];
%combine ex and in
p1=[ei_cp ;ei_ct ;ei_cci];p2=[ei_pv ;ei_som ;ei_vip];
par=[];par=[p1' p2']';
g1=[];g1=ones(1,length(p1)); 
g2=[];g2=ones(1,length(p2))*2;
gro=[];gro=[g1 g2]';
fig7= figure;set(fig7, 'Name', 'Barplot groups');set(fig7, 'Position', [400, 500, 250, 250]);set(gcf,'color','w');
hold on;line([0 3],[0 0],'LineStyle',':','Color','k','LineWidth',1)
violins = violinplot(par, gro,'ViolinColor',[[0.3 0.3 0.3]; [0.6 0.6 0.6]],'ShowMean', true);box off;
xlim([0 3]);ylabel('EX-IN index');set(gca,'FontSize',12);
% h = gca;h.XAxis.Visible = 'off';
% set(gca,'xtick',[]);
xticklabels({'PNs','INs'});
ylim([-1.2 1.2]);
%statitics
p_ei=ranksum([ei_cp ;ei_ct ;ei_cci],[ei_pv ;ei_som ;ei_vip]);
set(gca,'FontSize',11);set(gca,'TickDir','out'); hold on;text(1.25,1,['**'],'FontSize',18);
%% example 
temp1=[];ntsr_cs1=[];
for i=1:6
temp1(i,:)=cell_selecter(Ephys,'drugs',0,'label',4,'geno',6,'sol',2,'qualityinput',1,'pair',i);
end
ntsr_cs1=sum(temp1);
temp1=[];ntsr_cs2=[];
for i=1:6
temp1(i,:)=cell_selecter(Ephys,'drugs',1,'label',4,'geno',6,'sol',2,'qualityinput',1,'pair',i);
end
ntsr_cs2=sum(temp1);
ntsr_cs=[];
ntsr_cs=ntsr_cs1+ntsr_cs2;
disp([num2str(sum(ntsr_cs)) ' ntsr_cs cells, nr ' num2str(find(ntsr_cs==1))])

%CS NTSR1 paired (no drugs and before washin) CORRESPONDING CPN cells
temp1=[];ntsr_cpncs1=[];
for i=1:6
temp1(i,:)=cell_selecter(Ephys,'drugs',0,'label',1,'geno',6,'sol',2,'qualityinput',1,'pair',i);
end
ntsr_cpncs1=sum(temp1);
%
temp1=[];ntsr_cpncs2=[];
for i=1:6
temp1(i,:)=cell_selecter(Ephys,'drugs',1,'label',1,'geno',6,'sol',2,'qualityinput',1,'pair',i);
end
ntsr_cpncs2=sum(temp1);
ntsr_cpncs=[];
ntsr_cpncs=ntsr_cpncs1+ntsr_cpncs2;
disp([num2str(sum(ntsr_cpncs)) ' ntsr_cpncs cells, nr ' num2str(find(ntsr_cpncs==1))])
%remove cell 253 from ntsr1 cause its pair is with TTX
ntsr_cs(253)=0;
disp([num2str(sum(ntsr_cs)) ' updated ntsr_cs cells, nr ' num2str(find(ntsr_cs==1))])
%% 
%CP
cnr=[];cnr=find(ntsr_cpncs==1)
a=4;
ex_trace=[];ex_trace=Ephys(cnr(a)).sub_traces_train(:,2);
in_trace=[];in_trace=Ephys(cnr(a)).sub_traces_train(:,1);
sample_rate=sr;
show_time=[3000:10000];
start_time=5000;
time_stop=2000;
yborders=[-500 1000];
tit='CP'

epsc_ipsc_example(ex_trace,in_trace,sample_rate,show_time,start_time,time_stop,yborders,tit);
%CT
cnr=[];cnr=find(ntsr_cs==1)
ex_trace=[];ex_trace=Ephys(cnr(a)).sub_traces_train(:,2);
in_trace=[];in_trace=Ephys(cnr(a)).sub_traces_train(:,1);
sample_rate=sr;
show_time=[3000:10000];
start_time=5000;
time_stop=2000;
yborders=[-500 1000];
tit='CT'

epsc_ipsc_example(ex_trace,in_trace,sample_rate,show_time,start_time,time_stop,yborders,tit);
%% 

cnr=1;
ov_min=-250;ov_max=800;
temp=[];
temp=find(ntsr_cpncs_idx==1);
fig4=figure;set(fig4, 'Position', [200, 800, 400, 200]);set(gcf,'color','w');
subplot(1,2,1)
plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,1),'Color','r','LineWidth',1);set(gca,'box','off');
hold on;plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,2),'Color','b','LineWidth',1);set(gca,'box','off');hold on;ylim([ov_min-10 ov_max]);
hold on;plot([0.25*sr 0.25*sr],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
%hold on;title('CP');
axis off;

cnr=3
subplot(1,2,2)
temp=[];
temp=find(ntsr_cs_idx==1);
plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,1),'Color','b','LineWidth',1);set(gca,'box','off');
hold on;plot(Ephys(temp(cnr)).sub_traces_train(1:1*sr,2),'Color','r','LineWidth',1);set(gca,'box','off');
hold on;ylim([ov_min-10 ov_max]);
%title('Ntsr1','Color','k');
hold on;plot([0.25*sr 0.25*sr],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
axis off;
%Scale bar
 scale_x= 200;
 scale_y= 100;
 %scale barx
 hold on;x1= (1250*srF)-(100*srF);x2=1250*srF;p1=plot([x1 x2],[ov_min-10 ov_min-10],'-','Color','k','LineWidth',1.5);
 %scale bary
 hold on;y2= (ov_min-10)+scale_y;y1=(ov_min-10);p2=plot([x2 x2],[y1 y2],'-','Color','k','LineWidth',1.5);
