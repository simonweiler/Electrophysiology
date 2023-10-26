function [data] = timing_readout(Ephys, cell_filter, e_i_ratio, fc_1, sub_d,stimuli_type,rsp_window)
%function to read out timing of EPSC and ISPC in the same cell
temp=[];t_ex=[];t_in=[];t_in_ex=[];trace_smooth_ex=[];trace_smooth_in=[];trace_smooth_in_ex=[];
temp=find(cell_filter==1);sr=20000;
close all;
%% Check whic trace is EPSC and which one is IPSC and then run 

for i=1:length(temp)
    stimuli_trace=[];
     if stimuli_type==1
        stimuli_trace=Ephys(temp(i)).sub_traces_train;
        elseif stimuli_type==2
        stimuli_trace=Ephys(temp(i)).sub_traces_high;
        else stimuli_type==3
        stimuli_trace=Ephys(temp(i)).sub_traces_highf;
        end
    try
       
    if max(stimuli_trace(1:1*sr,2))>max(stimuli_trace(1:1*sr,1))==1
[t_ex(i) trace_smooth_ex(i,:)]=time_to_peak_callosal(stimuli_trace(:,1),[3000:1:5000],rsp_window,20,0,fc_1,temp(i));
saveas(gcf, [ 'Cell EX ' num2str(temp(i)) '.jpg']);close all;
[t_in(i) trace_smooth_in(i,:)]=time_to_peak_callosal(stimuli_trace(:,2),[3000:1:5000],rsp_window,20,1,fc_1,temp(i));
saveas(gcf, [ 'Cell IN ' num2str(temp(i)) '.jpg']);close all;
[t_in_ex(i) trace_smooth_in_ex(i,:)]=time_to_peak_callosal(stimuli_trace(:,2),[3000:1:5000],rsp_window,20,0,fc_1,temp(i));
saveas(gcf, [ 'Cell exin' num2str(temp(i)) '.jpg']);close all;
cell_id_save(i)=temp(i);
    else
      [t_ex(i) trace_smooth_ex(i,:)]=time_to_peak_callosal(stimuli_trace(:,2),[3000:1:5000],rsp_window,20,0,fc_1,temp(i));
      saveas(gcf, [ 'Cell EX ' num2str(temp(i)) '.jpg']);close all;
      [t_in(i) trace_smooth_in(i,:)]=time_to_peak_callosal(stimuli_trace(:,1),[3000:1:5000],rsp_window,20,1,fc_1,temp(i));
      saveas(gcf, [ 'Cell IN ' num2str(temp(i)) '.jpg']);close all;
      [t_in_ex(i) trace_smooth_in_ex(i,:)]=time_to_peak_callosal(stimuli_trace(:,1),[3000:1:5000],rsp_window,20,0,fc_1,temp(i));
      saveas(gcf, [ 'Cell exin ' num2str(temp(i)) '.jpg']);close all;
      cell_id_save(i)=temp(i);
    end
    catch
     t_ex(i)=NaN;
     t_in(i)=NaN; 
     t_in_ex(i)=NaN;   
     cell_id_save(i)=temp(i);
    end
   
end

%% PAIRED COMPARISON
%delete out nonsense 
%% Paired comparison of onset latency between EPSC and IPSC  panel d
if sub_d==1
 t_ex_sub=[];t_ex_sub=t_ex(find(e_i_ratio>0));
 t_in_sub=[];t_in_sub=t_in(find(e_i_ratio>0));
 t_in_ex_sub=[];t_in_ex_sub=t_in_ex(find(e_i_ratio>0));
 cell_id_sub=[];cell_id_sub=cell_id_save(find(e_i_ratio>0));
 
 subex=[];subex=find(t_in_ex_sub<t_ex_sub);
 t_ex_sub(subex)=t_in_ex_sub(subex); 

t_ein=[];
 
t_ein=[t_ex_sub' t_in_sub' cell_id_sub'];

else
  t_ein=[];
  t_ex_sub=t_ex;
   t_in_sub=t_in;
 t_ein=[t_ex_sub' t_in_sub'];
end
% [gh gm]=find(isnan(t_ein));
% t_ein(unique(gh),:)=[];


data=[];data=t_ein;

end
