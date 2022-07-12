function timing_readout(Ephys, cell_filter, e_i_ratio, fc_1, sub_d)

temp=[];t_ex=[];t_in=[];t_in_ex=[];trace_smooth_ex=[];trace_smooth_in=[];trace_smooth_in_ex=[];

temp=find(cell_filter==1);
sr=20000;
close all
for i=1:length(temp)
    try
    if max(Ephys(temp(i)).sub_traces_train(1:1*sr,2))>max(Ephys(temp(i)).sub_traces_train(1:1*sr,1))==1
[t_ex(i) trace_smooth_ex(i,:)]=time_to_peak(Ephys(temp(i)).sub_traces_train(:,1),[3000:1:5000],[5000:1:6000],20,0,fc_1);
[t_in(i) trace_smooth_in(i,:)]=time_to_peak(Ephys(temp(i)).sub_traces_train(:,2),[3000:1:5000],[5000:1:6000],20,1,fc_1);
[t_in_ex(i) trace_smooth_in_ex(i,:)]=time_to_peak(Ephys(temp(i)).sub_traces_train(:,2),[3000:1:5000],[5000:1:6000],20,0,fc_1);
    else
      [t_ex(i) trace_smooth_ex(i,:)]=time_to_peak(Ephys(temp(i)).sub_traces_train(:,2),[3000:1:5000],[5000:1:6000],20,0,fc_1);
      [t_in(i) trace_smooth_in(i,:)]=time_to_peak(Ephys(temp(i)).sub_traces_train(:,1),[3000:1:5000],[5000:1:6000],20,1,fc_1);
      [t_in_ex(i) trace_smooth_in_ex(i,:)]=time_to_peak(Ephys(temp(i)).sub_traces_train(:,1),[3000:1:5000],[5000:1:6000],20,0,fc_1);
    end
    catch
     t_ex(i)=NaN;
     t_in(i)=NaN; 
     t_in_ex(i)=NaN;   
    end
end
close all;
%% 

close all;%% Paired compariosn of onset latency between EPSC and IPSC  panel d
if sub_d==1
    t_ex_sub=[];t_ex_sub=t_ex(find(e_i_ratio>0));
 t_in_sub=[];t_in_sub=t_in(find(e_i_ratio>0));
 t_in_ex_sub=[];t_in_ex_sub=t_in_ex(find(e_i_ratio>0));
 subex=[];subex=find(t_in_ex_sub<t_ex_sub);
 t_ex_sub(subex)=t_in_ex_sub(subex);
t_ein=[];
t_ein=[t_ex_sub' t_in_sub'];
else
  t_ein=[];
  t_ex_sub=t_ex;
   t_in_sub=t_in;
 t_ein=[t_ex_sub' t_in_sub'];
end
[gh gm]=find(isnan(t_ein));
t_ein(unique(gh),:)=[];
cl={'r','b'};
data=[];data=t_ein;;
paired_plot_box(data,cl);
xticklabels({'EX','IN'});ylabel('Onset Latency (ms)');set(gca,'FontSize',10);
%statistics 
%test for normality: 
kstest(t_ex_sub');
kstest(t_in_sub');
%pvalue of paired signrank test:  4.8828e-04
[p1]=signrank(t_ein(:,1),t_ein(:,2));
end
