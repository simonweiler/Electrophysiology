%% Plot different light evoked EPSC/EPSP/IPSC of group of cells 

%% DEPENDENCIES
    %uipickfiles.m
    %cell_selecter.m
    %use plot_slice_opto.m
%% LEGENDS
%genotypes 1=Tlx3Cre; 2=Rpb4Cre; 3=Syt6Cre; 4=GADcre tdtomato; 5= PVcre
%tdtomat; 6=Ntsr1cre tdtomato ;7=Penk; 8=SOM; 9=VIP;

%labels
%1=callosal; 2=GAD; 3=PV; 4=Ntsr1 tdtomato; 5= Penk
%6=SOM; VIP=8

% Solution 
%K-gluc=1
%Cs-gluc=2
%% Paths, modify accordingly
CP_K_gluc='D:\Postdoc_Margrie\Projects\Callosal\output\CP_K_gluc';
CT_K_gluc='D:\Postdoc_Margrie\Projects\Callosal\output\CT_K_gluc';
GAD_K_gluc='D:\Postdoc_Margrie\Projects\Callosal\output\GAD_K_gluc';
PV_K_gluc='D:\Postdoc_Margrie\Projects\Callosal\output\PV_K_gluc';
SOM_K_gluc='D:\Postdoc_Margrie\Projects\Callosal\output\SOM_K_gluc';
VIP_K_gluc='D:\Postdoc_Margrie\Projects\Callosal\output\VIP_K_gluc';
PENK_K_gluc='D:\Postdoc_Margrie\Projects\Callosal\output\PENK_K_gluc';
UN_K_gluc='D:\Postdoc_Margrie\Projects\Callosal\output\UN_K_gluc';
%% Load data structure 
str_L6    = 'D:\Postdoc_Margrie\Projects\Callosal\output';
folder_list = uipickfiles('FilterSpec',str_L6);
load(char(folder_list));
%sampling rate
srF=20;sr=20000;
%% use plot_slice_opto.m 
%% 1) CP K-gluc all (111 cells)
cp_k=cell_selecter(Ephys,'drugs',0,'label',1,'sol',1);
cd(CP_K_gluc);
temp=[];
temp=find(cp_k==1);
for i=1:length(temp)
    plot_slice_opto(Ephys, temp(i), 1,sr);
    saveas(gcf, ['Cell ' num2str(temp(i)) '.jpg']);close all;
end
%% 2) CT K-gluc all (26 cells)
temp1=[];temp1=cell_selecter(Ephys,'drugs',0,'label',4,'sol',1,'geno',6);
temp2=[];temp2=cell_selecter(Ephys,'drugs',1,'label',4,'sol',1,'geno',6);
ct_k=temp1+temp2;
cd(CT_K_gluc);
temp=[];
temp=find(ct_k==1);
for i=1:length(temp)
    plot_slice_opto(Ephys, temp(i), 1,sr);
    saveas(gcf, ['Cell ' num2str(temp(i)) '.jpg']);close all;
end
%% 3) GAD K-gluc all (30 cells)
gad_k=cell_selecter(Ephys,'drugs',0,'label',2,'sol',1,'geno',4);
cd(GAD_K_gluc);
temp=[];
temp=find(gad_k==1);
for i=1:length(temp)
    plot_slice_opto(Ephys, temp(i), 1,sr);
    saveas(gcf, ['Cell ' num2str(temp(i)) '.jpg']);close all;
end
%% 4) PV K-gluc all (14)
pv_k=cell_selecter(Ephys,'drugs',0,'label',3,'sol',1,'geno',5);
cd(PV_K_gluc);
temp=[];
temp=find(pv_k==1);
for i=1:length(temp)
    plot_slice_opto(Ephys, temp(i), 1,sr);
    saveas(gcf, ['Cell ' num2str(temp(i)) '.jpg']);close all;
end
%% 5) SST K-gluc all (33)
som_k=cell_selecter(Ephys,'drugs',0,'label',6,'sol',1,'geno',8);
cd(SOM_K_gluc);
temp=[];
temp=find(som_k==1);
for i=1:length(temp)
    plot_slice_opto(Ephys, temp(i), 1,sr);
    saveas(gcf, ['Cell ' num2str(temp(i)) '.jpg']);close all;
end
%% 6) VIP K-gluc all (8)
vip_k=cell_selecter(Ephys,'drugs',0,'label',8,'sol',1,'geno',9);
cd(VIP_K_gluc);
temp=[];
temp=find(vip_k==1);
for i=1:length(temp)
    plot_slice_opto(Ephys, temp(i), 1,sr);
    saveas(gcf, ['Cell ' num2str(temp(i)) '.jpg']);close all;
end
%% 7) PENK K
penk_k=cell_selecter(Ephys,'drugs',0,'label',5,'sol',1,'geno',7);
cd(PENK_K_gluc);
temp=[];
temp=find(penk_k==1);
for i=1:length(temp)
    plot_slice_opto(Ephys, temp(i), 1,sr);
    saveas(gcf, ['Cell ' num2str(temp(i)) '.jpg']);close all;
end
%% 8) unlabelled K
un_k=cell_selecter(Ephys,'drugs',0,'label',0,'sol',1);
cd(UN_K_gluc);
temp=[];
temp=find(un_k==1);
for i=1:length(temp)
    plot_slice_opto(Ephys, temp(i), 1,sr);
    saveas(gcf, ['Cell ' num2str(temp(i)) '.jpg']);close all;
end
%% CS solution 
CP_CS_gluc='D:\Postdoc_Margrie\Projects\Callosal\output\CP_CS_gluc';
CT_CS_gluc='D:\Postdoc_Margrie\Projects\Callosal\output\CT_CS_gluc';
GAD_CS_gluc='D:\Postdoc_Margrie\Projects\Callosal\output\GAD_CS_gluc';
PV_CS_gluc='D:\Postdoc_Margrie\Projects\Callosal\output\PV_CS_gluc';
SOM_CS_gluc='D:\Postdoc_Margrie\Projects\Callosal\output\SOM_CS_gluc';
VIP_CS_gluc='D:\Postdoc_Margrie\Projects\Callosal\output\VIP_CS_gluc';
PENK_CS_gluc='D:\Postdoc_Margrie\Projects\Callosal\output\PENK_CS_gluc';
%% CP_CS (77 cells)
temp1=[];temp1=cell_selecter(Ephys,'drugs',0,'label',1,'sol',2);
temp2=[];temp2=cell_selecter(Ephys,'drugs',1,'label',1,'sol',2);
cp_cs=temp1+temp2;
cd(CP_CS_gluc);
temp=[];
temp=find(cp_cs==1);
for i=1:length(temp)
    plot_slice_opto(Ephys, temp(i), 2,sr);
    saveas(gcf, ['Cell ' num2str(temp(i)) '.jpg']);close all;
end
%% CT_CS (20)
temp1=[];temp1=cell_selecter(Ephys,'drugs',0,'label',4,'sol',2,'geno',6);
temp2=[];temp2=cell_selecter(Ephys,'drugs',1,'label',4,'sol',2,'geno',6);
ct_cs=temp1+temp2;
cd(CT_CS_gluc);
temp=[];
temp=find(ct_cs==1);
for i=1:length(temp)
    plot_slice_opto(Ephys, temp(i), 2,sr);
    saveas(gcf, ['Cell ' num2str(temp(i)) '.jpg']);close all;
end
%% GAD_CS (1)
gad_cs=cell_selecter(Ephys,'drugs',0,'label',2,'sol',2,'geno',4);
cd(GAD_CS_gluc);
temp=[];
temp=find(gad_cs==1);
for i=1:length(temp)
    plot_slice_opto(Ephys, temp(i), 1,sr);
    saveas(gcf, ['Cell ' num2str(temp(i)) '.jpg']);close all;
end
%% PV_CS
pv_cs=cell_selecter(Ephys,'drugs',0,'label',3,'sol',2,'geno',5);
cd(PV_CS_gluc);
temp=[];
temp=find(pv_cs==1);
for i=1:length(temp)
    plot_slice_opto(Ephys, temp(i), 1,sr);
    saveas(gcf, ['Cell ' num2str(temp(i)) '.jpg']);close all;
end
%% %% SST_CS
som_cs=cell_selecter(Ephys,'drugs',0,'label',6,'sol',2,'geno',8);
cd(SOM_CS_gluc);
temp=[];
temp=find(som_cs==1);
for i=1:length(temp)
    plot_slice_opto(Ephys, temp(i), 1,sr);
    saveas(gcf, ['Cell ' num2str(temp(i)) '.jpg']);close all;
end
%% VIP_CS
vip_cs=cell_selecter(Ephys,'drugs',0,'label',8,'sol',2,'geno',9);
cd(VIP_CS_gluc);
temp=[];
temp=find(vip_cs==1);
for i=1:length(temp)
    plot_slice_opto(Ephys, temp(i), 1,sr);
    saveas(gcf, ['Cell ' num2str(temp(i)) '.jpg']);close all;
end
%% PENK_CS
penk_cs=cell_selecter(Ephys,'drugs',0,'label',5,'sol',2,'geno',7);
cd(PENK_CS_gluc);
temp=[];
temp=find(penk_cs==1);
for i=1:length(temp)
    plot_slice_opto(Ephys, temp(i), 1,sr);
    saveas(gcf, ['Cell ' num2str(temp(i)) '.jpg']);close all;
end
