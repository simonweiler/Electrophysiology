
%% Load data structure 
str_L6    = 'D:\Postdoc_Margrie\Projects\Callosal\output';
folder_list = uipickfiles('FilterSpec',str_L6);
load(char(folder_list));
%sampling rate
srF=20;
sr=20000;
% 
%% 

%K
temp1=[];som_k=[];
temp1=cell_selecter(Ephys,'drugs',0,'label',6,'geno',8,'sol',1,'optovariant',1);
som_k=temp1;
[epsp_som] = readout_amp_epsp(Ephys,som_k,1,sr);
%% 

figure;
   temp=find(som_k==1);
   hold on;
for i=1:length(find(som_k==1)) 
    try
plot(ones(1,length(temp)),max(Ephys(temp(i)).IV.spikecount),'*k');
    catch
    end
end
hold on;
temp2=find(epsp_som>30);
for i=1:length(find(epsp_som>30)) 
    try
plot(ones(1,length(temp2))*2,max(Ephys(temp(temp2(i))).IV.spikecount),'*r');
m(i)=max(Ephys(temp(temp2(i))).IV.spikecount);
    catch
    end
end
hold on;
temp3=find(epsp_som<30);
for i=1:length(find(epsp_som<30)) 
    try
plot(ones(1,length(temp3))*3,max(Ephys(temp(temp3(i))).IV.spikecount),'*b');
    catch
    end
end
%% 
figure;
   temp=find(som_k==1);
   hold on;
for i=1:length(find(som_k==1)) 
    try
plot(ones(1,length(temp)),max(Ephys(temp(i)).Passive.tau),'*k');
    catch
    end
end
hold on;
temp2=find(epsp_som>30);
for i=1:length(find(epsp_som>30)) 
    try
plot(ones(1,length(temp2))*2,max(Ephys(temp(temp2(i))).Passive.tau),'*r');
m(i)=max(Ephys(temp(temp2(i))).IV.spikecount);
    catch
    end
end
hold on;
temp3=find(epsp_som<30);
for i=1:length(find(epsp_som<30)) 
    try
plot(ones(1,length(temp3))*3,max(Ephys(temp(temp3(i))).Passive.tau),'*b');
    catch
    end
end
