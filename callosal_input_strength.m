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
%% readout input strengths
[ntsr_cs ntsr_cpncs ntsr_k ntsr_cpnk ...
    penk_cs penk_cpncs penk_k penk_cpnk ...
    som_cs som_cpncs som_k som_cpnk ...
    pv_cs pv_cpncs pv_k pv_cpnk] = readout_input_strength_callosal(Ephys);
%% cs and k ratio long pulse 
ntsr_ratios=[];
for i=1:size(ntsr_cs,2)
ntsr_er=[];
ntsr_er=[ntsr_cs(:,i)./ntsr_cpncs(:,i) ;ntsr_k(:,i)./ntsr_cpnk(:,i)];
ntsr_er(find(ntsr_er==Inf))=NaN;
ntsr_ratios(:,i)=ntsr_er;
end

penk_ratios=[];
for i=1:size(penk_cs,2)
penk_er=[];
penk_er=[penk_cs(:,i)./penk_cpncs(:,i) ;penk_k(:,i)./penk_cpnk(:,i)];
penk_er(find(penk_er==Inf))=NaN;
penk_ratios(:,i)=penk_er;
end

som_ratios=[];
for i=1:size(som_cs,2)
som_er=[];
som_er=[som_cs(:,i)./som_cpncs(:,i) ;som_k(:,i)./som_cpnk(:,i)];
som_er(find(som_er==Inf))=NaN;
som_ratios(:,i)=som_er;
end

pv_ratios=[];
for i=1:size(pv_cs,2)
pv_er=[];
pv_er=[pv_cs(:,i)./pv_cpncs(:,i) ;pv_k(:,i)./pv_cpnk(:,i)];
pv_er(find(pv_er==Inf))=NaN;
pv_ratios(:,i)=pv_er;
end
%% Barplots for all 
fig1= figure;set(fig1, 'Name', 'Barplot groups');set(fig1, 'Position', [400, 500, 500, 200]);set(gcf,'color','w');
data=[];data=ntsr_ratios;
for i=1:size(ntsr_cs,2)
hold on;
plot(ones(1,length(data(:,i)))*i,data(:,i),'ko','MarkerEdgeColor',[0.7,0.7,0.7],'MarkerSize',3);
hold on;
b2=bar(i,nanmean(data(:,i)));b2.FaceColor='r';
hold on;
er=errorbar(i,nanmean(data(:,i)),nanstd(data(:,i))/sqrt(size(data(:,i),1)));er.Color = [0 0 0];er.LineWidth=1.5;er.LineStyle = 'none'; hold on;
end
xlim([0 10]);
plot(0:i,ones(1,size(data,2)+1))

