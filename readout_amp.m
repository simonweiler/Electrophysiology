function [epsc ipsc e_i_ratio] = readout_amp(str,cells_idx,stim_type,sol,ind1,ind2)
%read out the overall maximum amplitude for EPSC and ISPC #SW211118
%% OUTPUT
%epsc=max epsc amplitude
%ipsc=max ipsc amplitude
%e_i_ratio= E/Iratio
%% INPUT
%str: Ephys structure
%cells_idx : cells id of desired cells
%stim_type : 1=train; 2=highfr 1;3=highfr2
sr=20000;
%% function 
temp=find(cells_idx==1);
%1:read out train stimulus (5 long pulses)
if stim_type==1   
    for i=1:length(find(cells_idx==1));  
     
       if isempty(str(temp(i)).train_n)==0;
           if sol==2
               try
          if max(str(temp(i)).sub_traces_train(1:1*sr,ind2))>max(str(temp(i)).sub_traces_train(1:1*sr,ind1))==1
         epsc(i)=max(abs(str(temp(i)).train_n(:,ind1)));
         ipsc(i)=max(abs(str(temp(i)).train_p(:,ind2)));
          else
         epsc(i)=max(abs(str(temp(i)).train_n(:,ind2)));
         ipsc(i)=max(abs(str(temp(i)).train_p(:,ind1)));
          end
              catch 
         epsc(i)=max(abs(str(temp(i)).train_n(:)));
         ipsc(i)=max(abs(str(temp(i)).train_p(:))); 
             end
          else sol==1
         epsc(i)=max(abs(str(temp(i)).train_n(:)));
         ipsc(i)=max(abs(str(temp(i)).train_p(:)));  
           end
       else
           epsc(i)=NaN;
           ipsc(i)=NaN;
       end
    end

%2:read out high freq stimulus (25, short pulses)   
elseif stim_type==2
    for i=1:length(find(cells_idx==1));   
      if isempty(str(temp(i)).high_n)==0;
            if sol==2
             try
     if max(str(temp(i)).sub_traces_high(1:1*sr,ind2))>max(str(temp(i)).sub_traces_high(1:1*sr,ind1))==1
         epsc(i)=max(abs(str(temp(i)).high_n(:,ind1)));
         ipsc(i)=max(abs(str(temp(i)).high_p(:,ind2)));
          else
         epsc(i)=max(abs(str(temp(i)).high_n(:,ind2)));
         ipsc(i)=max(abs(str(temp(i)).high_p(:,ind1)));
     end
             catch 
         epsc(i)=max(abs(str(temp(i)).high_n(:)));
         ipsc(i)=max(abs(str(temp(i)).high_p(:))); 
             end
            else sol==1
         epsc(i)=max(abs(str(temp(i)).high_n(:)));
         ipsc(i)=max(abs(str(temp(i)).high_p(:)));  
            end
    
        else
         epsc(i)=NaN;
         ipsc(i)=NaN;  
        end
    end
%3:read out highest freq stimulus (50, short pulses)   
else stim_type==3
    for i=1:length(find(cells_idx==1)); 
        if isempty(str(temp(i)).highf_n)==0;
        if sol==2
            try
   if max(str(temp(i)).sub_traces_highf(1:1*sr,ind2))>max(str(temp(i)).sub_traces_highf(1:1*sr,ind1))==1
         epsc(i)=max(abs(str(temp(i)).highf_n(:,ind1)));
         ipsc(i)=max(abs(str(temp(i)).highf_p(:,ind2)));
          else
         epsc(i)=max(abs(str(temp(i)).highf_n(:,ind2)));
         ipsc(i)=max(abs(str(temp(i)).highf_p(:,ind1)));
   end
       catch 
         epsc(i)=max(abs(str(temp(i)).highf_n(:)));
         ipsc(i)=max(abs(str(temp(i)).highf_p(:))); 
             end
        else sol==1
            epsc(i)=max(abs(str(temp(i)).highf_n(:)));
         ipsc(i)=max(abs(str(temp(i)).highf_p(:)));

        end
        else
         epsc(i)=NaN;
         ipsc(i)=NaN;  
        end
    end 
    
end
%% E / I ratio
for i=1:length(epsc);
  if epsc(i)==0 & ipsc(i)==0;
e_i_ratio(i)=NaN;
  else
      e_i_ratio(i)=epsc(i)/ipsc(i);
  end
 
end
e_i_ratio(find(isinf(e_i_ratio)))=NaN;
end