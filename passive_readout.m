function [rmp maxsp rheo rin tau sag trace spike_time] = passive_readout(data,idx)
%input: Ephys structure from ephys analysis 
%output: passive properties

temp=[];
for i=1:length(find(idx==1));
    temp=find(idx==1);
    if isempty(data(temp(i)).IV)==0 
    rmp(i)=data(temp(i)).IV.RMP; 
    maxsp(i)=max(data(temp(i)).IV.spikecount);
    else
        rmp(i)=NaN;
        maxsp(i)=NaN;
    end

 if isempty(data(temp(i)).Rheobase)==0 
    rheo(i)=data(temp(i)).Rheobase.rheo;
 elseif isempty(data(temp(i)).IV)==0 

      tr=[];tr=find(data(temp(i)).IV.spikecount>0);
      try
      rheo(i)=data(temp(i)).IV.stimvec(tr(1)); 
      catch
      rheo(i)=NaN;
      end
 else
      rheo(i)=NaN;
 end

     
  if  isempty(data(temp(i)).Passive)==0
    rin(i)=data(temp(i)).Passive.Rin;
    tau(i)=data(temp(i)).Passive.tau(1);
  else
    rin(i)=NaN;
    tau(i)=NaN;
  end

  if isempty(data(temp(i)).Sag)==0
   sag(i)=data(temp(i)).Sag.sagratio(1);
  else
   sag(i)=NaN;
  end


    md=[];
    if isempty(data(temp(i)).Rheobase)==0
     if isempty(find(data(temp(i)).Rheobase.spikecount(1:size(data(temp(i)).Rheobase.traces,2))>=2))==0
    md=find(data(temp(i)).Rheobase.spikecount(1:size(data(temp(i)).Rheobase.traces,2))>=2);
    trace(:,i)=data(temp(i)).Rheobase.traces(:,md(1));
     sp_time=[];
    diff_idx=[];
    diff_tr=[];
   diff_idx=find(diff(find(diff(trace(:,i))>0.5))>1);
   d_tr=find(diff(trace(:,i))>0.5);
   dif_tr=diff(find(diff(trace(:,i))>0.5));
    sp_time=d_tr(diff_idx+1);
    spike_time(i)=sp_time(1);
   else
      md=find(data(temp(i)).Rheobase.spikecount(1:size(data(temp(i)).Rheobase.traces,2))>=1);
    trace(:,i)=data(temp(i)).Rheobase.traces(:,md(1));
    sp_time=[];
    diff_idx=[];
    diff_tr=[];
   diff_idx=find(diff(find(diff(trace(:,i))>0.5))>1);
   d_tr=find(diff(trace(:,i))>0.5);
   dif_tr=diff(find(diff(trace(:,i))>0.5));
    sp_time=d_tr(1);
    spike_time(i)=sp_time(1);
     end
      
    elseif isempty(data(temp(i)).IV)==0
         tr=[];tr=find(data(temp(i)).IV.spikecount>0);
         try
         trace(:,i)=data(temp(i)).IV.traces(:,tr(1)); 
         catch
         trace(:,i)=ones(25000,1)*NaN;   
         end

    else
        trace(:,i)=ones(25000,1)*NaN;
    end
end

end