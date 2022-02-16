function [epsp] = readout_amp_epsp(str,cells_idx,stim_type,sr)
%read out the overall maximum amplitude for EPSC and ISPC #SW211118
%% OUTPUT
%epsc=max epsc amplitude
%ipsc=max ipsc amplitude
%e_i_ratio= E/Iratio
%% INPUT
%str: Ephys structure
%cells_idx : cells id of desired cells
%stim_type : 1=train; 2=highfr 1;3=highfr2
%% function 
temp=find(cells_idx==1);
%1:read out train stimulus (5 long pulses)
if stim_type==1   
    for i=1:length(find(cells_idx==1));  
       if isempty(str(temp(i)).train_p)==0;
           amp=[];
           for k=1:size(str(temp(i)).sub_traces_train,2)
           if isempty(find(str(temp(i)).sub_traces_train(6.24*sr:end,k)<-50))==0
               amp(k)=NaN;
           else isempty(find(str(temp(i)).sub_traces_train(6.24*sr:end,k)<-50))==1
               amp(k)=max(abs(str(temp(i)).train_p(:,k)));
           end
           end
   epsp(i)=max(amp(find(~isnan(amp))));
       else
           epsp(i)=NaN;
          
       end
    end
%2:read out high freq stimulus (25, short pulses)   
elseif stim_type==2
     for i=1:length(find(cells_idx==1));  
       if isempty(str(temp(i)).high_p)==0;
           amp=[];
           for k=1:size(str(temp(i)).sub_traces_high,2)
           if isempty(find(str(temp(i)).sub_traces_high(6.24*sr:end,k)<-50))==0
               amp(k)=NaN;
           else isempty(find(str(temp(i)).sub_traces_high(6.24*sr:end,k)<-50))==1
               amp(k)=max(abs(str(temp(i)).high_p(:,k)));
           end
           end
   epsp(i)=max(amp(find(~isnan(amp))));
       else
           epsp(i)=NaN;
          
       end
    end
%3:read out highest freq stimulus (50, short pulses)   
else stim_type==3
      for i=1:length(find(cells_idx==1));  
       if isempty(str(temp(i)).highf_p)==0;
           amp=[];
           for k=1:size(str(temp(i)).sub_traces_highf,2)
           if isempty(find(str(temp(i)).sub_traces_highf(6.24*sr:end,k)<-50))==0
               amp(k)=NaN;
           else isempty(find(str(temp(i)).sub_traces_highf(6.24*sr:end,k)<-50))==1
               amp(k)=max(abs(str(temp(i)).highf_p(:,k)));
           end
           end
   epsp(i)=max(amp(find(~isnan(amp))));
       else
           epsp(i)=NaN;
          
       end
end

end