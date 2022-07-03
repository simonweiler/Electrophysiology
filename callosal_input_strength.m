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
thr_ratio=15;
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
pv_er(find(pv_er==Inf))=NaN;
pv_ratios(:,i)=pv_er;
end
%outlier removal
idr=[];idc=[];
[idr idc] = find(pv_ratios>thr_ratio);
for i=1:length(idc)
pv_ratios(idr(i),idc(i))=NaN;
end
%% Barplots for all 
typ_psc=1;
cl={'m','r','c','b'};

data={};data={ntsr_ratios(:,typ_psc) penk_ratios(:,typ_psc) som_ratios(:,typ_psc) pv_ratios(:,typ_psc)};
for i=1:length(data)
categ{:,i}=[ones(1,length(data{:,i}))*i]';
end
con_data=vertcat(data{:});
con_categ=vertcat(categ{:});
%figure
fig1= figure;set(fig1, 'Name', 'Barplot groups');set(fig1, 'Position', [400, 500, 150, 180]);set(gcf,'color','w');
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
plot(0:i+1,ones(1,size(data,2)+2),'LineStyle',':','Color','k');
xticks(1:4);
xticklabels({'Ntsr','Penk','SST','PV'});