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
    stimtcur=[];stimtcur=[Ephys(i).pulse_specs.amp_current];
    for k=1:length(stimty)
   
    traces_during=[Ephys(i).spikelight_traces(:,1:5,find(stimty==4 & stimtcur>40))];  
        
      end
   traces_all{:,i}= traces_during;

    end



%stim=horzcat(stimtype{:})';
%% just to browse 
m=3;cellnr=2;
trace_temp1=[];
trace_temp2=[];  
trace_temp1=traces_all{1,cellnr} (:,1,m);
trace_temp2=traces_all{1, cellnr} (:,1,m+1);
figure;plot(trace_temp1(:,1),'b');
hold on;plot(trace_temp2(:,1),'k');

%% 


%% for on and off stimuli
m=1;
trace_temp1=[];
trace_temp2=[];  
trace_temp1=traces_all{1,12} (:,1:5,m);
trace_temp2=traces_all{1, 12} (:,1:5,m+1);
%figure;plot(trace_temp1(15000:18000,1:5),'b');hold on;plot(trace_temp2(15000:18000,1:5),'k');
figure;plot(trace_temp1(:,1),'b');hold on;plot(trace_temp2(:,1),'k');
%% 
m=1;
for i=1:length(Ephys)
  tr1=[];tr2=[];
  trace_temp1=[];
trace_temp2=[];
  try
trace_temp1=traces_all{1, i} (15000:18000,1:5,m);
trace_temp2=traces_all{1, i} (15000:18000,1:5,m+1);


for k=1:5
[tr1] = findpeaks(trace_temp1(:,k),'MinPeakHeight',60);
if isempty(tr1)==1
count_tr1(k)=0;
else
 count_tr1(k)=length(tr1);   
end

[tr2] = findpeaks(trace_temp2(:,k),'MinPeakHeight',60);
if isempty(tr2)==1
count_tr2(k)=0;
else
 count_tr2(k)=length(tr2);   
end

end

spike_count1(i,:)=count_tr1;
spike_count2(i,:)=count_tr2;
spke_delta(i)=nanmean(count_tr1/0.15-count_tr2/0.15);

  catch
  end
end
