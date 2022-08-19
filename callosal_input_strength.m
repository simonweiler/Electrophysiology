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
%% readout input strengths for different groups
[ntsr_cs ntsr_cpncs ntsr_k ntsr_cpnk ...
    penk_cs penk_cpncs penk_k penk_cpnk ...
    som_cs som_cpncs som_k som_cpnk ...
    pv_cs pv_cpncs pv_k pv_cpnk] = readout_input_strength_callosal(Ephys);
%% cs and k ratio long pulse 
thr_ratio=10;
%NTSR
ntsr_ratios=[];
for i=1:size(ntsr_cs,2)
ntsr_er=[];
ntsr_er=[ntsr_cs(:,i)./ntsr_cpncs(:,i) ;ntsr_k(:,i)./ntsr_cpnk(:,i)];
ntsr_er(find(ntsr_er==Inf))=NaN;
ntsr_ratios(:,i)=ntsr_er;
end
%outlier removal
idr=[];idc=[];
[idr idc] = find(ntsr_ratios>thr_ratio);
for i=1:length(idc)
ntsr_ratios(idr(i),idc(i))=NaN;
end

%PENK
penk_ratios=[];
for i=1:size(penk_cs,2)
penk_er=[];
penk_er=[penk_cs(:,i)./penk_cpncs(:,i) ;penk_k(:,i)./penk_cpnk(:,i)];
penk_er(find(penk_er==Inf))=NaN;
penk_ratios(:,i)=penk_er;
end
%outlier removal
idr=[];idc=[];
[idr idc] = find(penk_ratios>thr_ratio);
for i=1:length(idc)
penk_ratios(idr(i),idc(i))=NaN;
end

%SOM
som_ratios=[];
for i=1:size(som_cs,2)
som_er=[];
som_er=[som_cs(:,i)./som_cpncs(:,i) ;som_k(:,i)./som_cpnk(:,i)];
som_er(find(som_er==Inf))=NaN;
som_ratios(:,i)=som_er;
end
%outlier removal
idr=[];idc=[];
[idr idc] = find(som_ratios>thr_ratio);
for i=1:length(idc)
som_ratios(idr(i),idc(i))=NaN;
end

%PV
pv_ratios=[];
for i=1:size(pv_cs,2)
pv_er=[];
pv_er=[pv_cs(:,i)./pv_cpncs(:,i) ;pv_k(:,i)./pv_cpnk(:,i)];
pv_er(pv_er==Inf)=NaN;
pv_ratios(:,i)=pv_er;
end
%outlier removal
idr=[];idc=[];
[idr idc] = find(pv_ratios>thr_ratio);
for i=1:length(idc)
pv_ratios(idr(i),idc(i))=NaN;
end
%% readout ntsr and ntsr cp index for plotting example 
temp1=[];ntsr_cs1=[];
for i=1:6
temp1(i,:)=cell_selecter(Ephys,'drugs',0,'label',4,'geno',6,'sol',2,'optovariant',1,'pair',i);
end
ntsr_cs1=sum(temp1);
temp1=[];ntsr_cs2=[];
for i=1:6
temp1(i,:)=cell_selecter(Ephys,'drugs',1,'label',4,'geno',6,'sol',2,'optovariant',1,'pair',i);
end
ntsr_cs2=sum(temp1);
ntsr_cs_idx=[];
ntsr_cs_idx=ntsr_cs1+ntsr_cs2;

%CS NTSR1 paired (no drugs and before washin) CORRESPONDING CPN cells
temp1=[];ntsr_cpncs1=[];
for i=1:6
temp1(i,:)=cell_selecter(Ephys,'drugs',0,'label',1,'geno',6,'sol',2,'optovariant',1,'pair',i);
end
ntsr_cpncs1=sum(temp1);
%
temp1=[];ntsr_cpncs2=[];
for i=1:6
temp1(i,:)=cell_selecter(Ephys,'drugs',1,'label',1,'geno',6,'sol',2,'optovariant',1,'pair',i);
end
ntsr_cpncs2=sum(temp1);
ntsr_cpncs_idx=[];
ntsr_cpncs_idx=ntsr_cpncs1+ntsr_cpncs2;
%remove cell 253 from ntsr1 cause its pair is with TTX
ntsr_cs_idx(253)=0
%% %% Plot exampe epsc / ipsc of retro and NTSR1
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

%% example paired box plot EPSC
tempd=[];tempd2=[];
col=1;
tempd=[ntsr_cpncs(:,col) ;ntsr_cpnk(:,col)];
tempd2=[ntsr_cs(:,col) ;ntsr_k(:,col)];
data=[];data=[tempd tempd2];
paired_plot_box(data,{'g','m'});
xticklabels({'CP','Ntsr1'});
ylabel('Light evoked EPSC amplitude');
%% example paired box plot IPSC
tempd=[];tempd2=[];
col=2;
tempd=[ntsr_cpncs(:,col) ;ntsr_cpnk(:,col)];
tempd2=[ntsr_cs(:,col) ;ntsr_k(:,col)];
data=[];data=[tempd tempd2];
paired_plot_box(data,{'g','m'});
xticklabels({'CP','Ntsr1'});
ylabel('Light evoked IPSC amplitude');
%% example paired box plot EI ratio
tempd=[];tempd2=[];
col=3;
tempd=[ntsr_cpncs(:,col) ;ntsr_cpnk(:,col)];
tempd2=[ntsr_cs(:,col) ;ntsr_k(:,col)];
data=[];data=[tempd tempd2];
paired_plot_box(data,{'g','m'});
xticklabels({'CP','Ntsr1'});
ylabel('E / I ratio');
ylim([0 1.2]);yticks(0:0.2:1.2);
 hold on;plot([0 3],[1 1],'--k');
%% Barplots for all 
typ_psc=3;
yl='E / I ratio';
cl={[0.7 0 0.4],[0.635 0.078 0.184],[0 0.5 1],[1 0.5 0]};

data={};data={ntsr_ratios(:,typ_psc) penk_ratios(:,typ_psc) som_ratios(:,typ_psc) pv_ratios(:,typ_psc)};
for i=1:length(data)
categ{:,i}=[ones(1,length(data{:,i}))*i]';
end
con_data=vertcat(data{:});
con_categ=vertcat(categ{:});
%figure
fig1= figure;set(fig1, 'Name', 'Barplot groups');set(fig1, 'Position', [400, 500, 230, 250]);set(gcf,'color','w');
for i=1:length(data)
hold on;
b2=bar(i,nanmean(data{:,i}));b2.FaceColor=cl{i};
hold on;
end
hold on;
plotSpread(data,'categoryIdx',con_categ,'categoryMarkers',{'.','.','.','.'},...
'categoryColors',{[0.7,0.7,0.7],[0.7,0.7,0.7],[0.7,0.7,0.7],[0.7,0.7,0.7]});
hold on;
for i=1:length(data)
er=errorbar(i,nanmean(data{:,i}),nanstd(data{:,i})/sqrt(size(data{:,i},1)));
er.Color = [0 0 0];er.LineWidth=1;er.LineStyle = 'none'; hold on;
end
xlim([0 5]);
ylim([0 10]);
plot(0:i+1,ones(1,size(data,2)+2),'LineStyle',':','Color','k','LineWidth',1.5);
xticks(1:4);
xticklabels({'Ntsr','Penk','SST','PV'});
ylabel({yl ; '(norm. by CP input)'});
set(gca,'FontSize',10);
%% Plot EI ratio histogram 
col=3;
data={};data={[ntsr_cpncs(:,col) ;penk_cpncs(:,col) ;som_cpncs(:,col) ;pv_cpncs(:,col)]...
    [ntsr_cs(:,col)] [penk_cs(:,col)] [som_cs(:,col)] [pv_cs(:,col)]};
 edges = [0:0.1:1.5];
  fig4=figure;set(fig4, 'Position', [200, 200, 240, 240]);set(gcf,'color','w');
  %for i=1:5
  %subplot(1,5,i);
hold on;h2=histogram(data{:,1},edges,'Normalization','probability');h2.FaceColor=[0 0.5 0.5];h2.EdgeColor='w';%h2.FaceAlpha=0.5;
 box off;xlabel({'E / I ratio' ; '(1 Hz stim freq.)'});ylabel('Relative counts');ylim([0 0.6]);
 hold on;plot([1 1],[0 0.5],'--k');
 hold on;plot([nanmedian(data{:,1}) nanmedian(data{:,1})],[0.5 0.5],...
     'Marker','v','MarkerFaceColor','k','MarkerEdgeColor','k');
 %title('25 Hz stim freq.','FontWeight','Normal');
;xlim([0 1.4]); set(gca,'FontSize',10)
 % end
  %% readout all cells type for timw windowing epsc ipsc
  temp2=[];cell_filter=[];
  rd=[0 1];
for i=1:2
%temp2(i,:)=cell_selecter(Ephys,'drugs',rd(i),'label',1,'sol',2,'optovariant',1);
temp2(i,:)=cell_selecter(Ephys,'drugs',rd(i),'label',4,'geno',6,'sol',2,'optovariant',1);

end
cell_filter=sum(temp2);
%% 

%% 
out_psc=[];
[out_psc] = readout_amp_update(Ephys,cell_filter ,1,2,1,2);
  %% Time to peak ex and in for retro cells using first pulse of long train 
%CPN
fc_1=3.5;
timing_readout(Ephys, cell_filter, out_psc(:,3), fc_1, 1);

%% %NTSR1
fc_1=2.5;
timing_readout(Ephys, cell_filter, out_psc(:,3), fc_1, 0);
%% Penk
 fc_1=3;
timing_readout(Ephys, cell_filter, out_psc(:,3), fc_1, 1);
%% SOM
 fc_1=2.5;
timing_readout(Ephys, cell_filter, out_psc(:,3), fc_1, 0);
%% PV
fc_1=3;
timing_readout(Ephys, cell_filter, out_psc(:,3), fc_1, 1);
%% FEEDFOWRD INHIBITION
pv_k_cells=cell_selecter(Ephys,'drugs',0,'label',3,'sol',1);
epsp_k=[];
[epsp_k] = readout_amp_epsp(Ephys,pv_k_cells,2,sr);
%% CP
cpn_k_cells=cell_selecter(Ephys,'drugs',0,'label',1,'sol',1);
[epsp_cpn]=[];
[epsp_cpn] = readout_amp_epsp(Ephys,cpn_k_cells,1,sr);
%% SOM
som_k_cells=cell_selecter(Ephys,'drugs',0,'label',6,'geno',8,'sol',1,'optovariant',1);
[epsp_som]=[];
[epsp_som] = readout_amp_epsp(Ephys,som_k_cells,2,sr);
%% read out all GAD
in_k_cells=cell_selecter(Ephys,'drugs',0,'label',2,'sol',1);

%% 
temp=[];
temp=find(in_k_cells==1);
for i=1:length(temp)
    if  isempty(Ephys(temp(i)).IV)==0

maxspF(i)=max(Ephys(temp(i)).IV.spikecount);
    else
        maxspF(i)=NaN;
    end

end
slow_spiker=find(maxspF<50);%
nonFS_k_cells=temp(slow_spiker);
nonPV_k_cells=zeros(1,length(Ephys));
nonPV_k_cells(nonFS_k_cells)=1;
%% 

[epsp_nonPV]=[];
[epsp_nonPV] = readout_amp_epsp(Ephys,nonPV_k_cells,1,sr);
%% example cells
cnr=91;%
ov_min=-5;ov_max=100;
range=4000:7000;
temp=[];temp=find(cpn_k_cells==1);
fig4=figure;set(fig4, 'Position', [200, 200, 130, 300]);set(gcf,'color','w');
subplot(3,1,1)
plot(Ephys(temp(cnr)).sub_traces_high(range,2),'Color','k','LineWidth',1.3);set(gca,'box','off');
hold on;plot([1000 1000],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
%end
ylim([ov_min-10 ov_max]);title('CP','Color','k');
axis off;

cnr=1;%
ov_min=-5;ov_max=100;
temp=[];temp=find(som_k_cells==1);
subplot(3,1,2)
plot(Ephys(temp(cnr)).sub_traces_high(range,2),'Color','#A2142F','LineWidth',1.3);set(gca,'box','off');
%hold on;plot([1000 1000],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
%end
ylim([ov_min-10 ov_max]);title('nPV IN (SST)','Color','#A2142F');
axis off;

cnr=12;%
ov_min=-5;ov_max=100;
temp=[];temp=find(pv_k_cells==1);
subplot(3,1,3)
plot(Ephys(temp(cnr)).sub_traces_high(range,2),'Color',[0.8500 0.3250 0.0980],'LineWidth',1.3);set(gca,'box','off');
%hold on;plot([1000 1000],[ov_max ov_max],'Marker','v','MarkerFaceColor','c','MarkerEdgeColor','c');
%end
ylim([ov_min-10 ov_max]);title('PV IN','Color',[0.8500 0.3250 0.0980]);
%axis off;
%% 
epsp_nFS=[epsp_nonPV epsp_som];
data=[];data=[epsp_cpn epsp_nFS epsp_k]';
fig6= figure;set(fig6, 'Name', 'compare fraction spiking');set(fig6, 'Position', [200, 300, 200, 300]);set(gcf,'color','w');
 hold on;scatter(ones(length(epsp_cpn),1),epsp_cpn,25,'o','MarkerEdgeColor','k','MarkerFaceColor','k');
 hold on;scatter(ones(length(epsp_nFS),1)*2,epsp_nFS,25,'o','MarkerEdgeColor','k','MarkerFaceColor','#A2142F');
 hold on;scatter(ones(length(epsp_k),1)*3, epsp_k,25,'o','MarkerEdgeColor','k','MarkerFaceColor',[0.8500 0.3250 0.0980]);
 %error bars and mean of subthreshold
 temp=[];temp=epsp_k;
e=errorbar(1,nanmean(epsp_cpn),nanstd(epsp_cpn)/sqrt(sum(~isnan(epsp_cpn))));
e.Color = 'k';e.CapSize = 10;
  hold on;plot([0.6 1.4],[nanmean(epsp_cpn) nanmean(epsp_cpn)],'Color','k');
e=errorbar(2,nanmean(epsp_nFS),nanstd(epsp_nFS)/sqrt(sum(~isnan(epsp_nFS))));
e.Color = '#A2142F';e.CapSize = 10;
  hold on;plot([1.6 2.4],[nanmean(epsp_nFS) nanmean(epsp_nFS)],'Color','#A2142F');
  e=errorbar(3,nanmean(temp(find(temp<50))),nanstd(temp(find(temp<50)))/sqrt(sum(~isnan(temp(find(temp<50))))));
e.Color = [0.8500 0.3250 0.0980];e.CapSize = 10;
 hold on;plot([2.6 3.4],[nanmean(temp(find(temp<50))) nanmean(temp(find(temp<50)))],'Color',[0.8500 0.3250 0.0980]);
  hold on;p=plot([0 4],[50 50],':k');
 xlim([0 4]);xticks([1:1:3]);ylabel('Light evoked EPSP amplitude (mV)');xticklabels({'CP','nPV IN','PV IN'});xtickangle(45);
 set(gca,'FontSize',10);
 %ylim([0 35])
 breakyaxis([35 55]);
