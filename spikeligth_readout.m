%% Load data structure 
str_L6    = 'D:\Postdoc_Margrie\Projects\Callosal\output';
folder_list = uipickfiles('FilterSpec',str_L6);
load(char(folder_list));
%sampling rate
srF=20;
sr=20000;
%% 


for i=1:length(Ephys)
    stimty=[];stimty=[Ephys(i).pulse_specs.stimtype];
    
    for k=1:length(stimty)
   
    traces_during=[Ephys(i).spikelight_traces(:,1:5,find(stimty==3))];  
        
      end
   traces_all{:,i}= traces_during;

    end



%stim=horzcat(stimtype{:})';
%% 

%% 
m=1;
trace_temp1=[];
trace_temp2=[];
trace_temp1=traces_all{1, 2} (:,1:5,m);
trace_temp2=traces_all{1, 2} (:,1:5,m+1);
figure;plot(trace_temp1(:,1:5),'b');hold on;plot(trace_temp2(:,1:5),'k');