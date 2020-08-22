function [F2xRheo]=plot_intrinsic(data,idx,siz,srF,col,ty)
%plot various ephys traces
%data=input structure
%idx=cell filter idx (logical)
%siz=size for subplot
%col=color for traces
%ty=what type of trace (e.g. IV, Rheobase, Passive,Sag, Ramp)
sr=20000;
if ty==1
step=3;
for i=1:length(find(idx==1));
temp=find(idx==1);
 if isempty(data(temp(i)).IV)==0
     min_tr(:,i)=min(min(data(temp(i)).IV.traces(:,1:step:end))); 
     max_tr(:,i)=max(max(data(temp(i)).IV.traces(:,1:step:end))); 
    else
     min_tr(:,i)=NaN;  
 end
end
ov_min=nanmin(min_tr);
ov_max=nanmax(max_tr);
%% Plotting IV
fig1=figure;set(fig1, 'Position', [200, -200, 800, 800]);set(gcf,'color','w');
for i=1:length(find(idx==1));
temp=find(idx==1);
hold on;subplot(siz,siz,i);
if isempty(data(temp(i)).IV)==0
%Plot only every second trace atm for display
plot(data(temp(i)).IV.traces(:,1:step:end),'Color',col,'LineWidth',1);set(gca,'box','off');axis off;

hold on;ylim([ov_min-10 ov_max]);
hold on;text(-550*srF,data(temp(i)).IV.RMP,[num2str(data(temp(i)).IV.RMP),'mV'],'FontSize',9);
 %Scale bar
 scale_x= 200;
 scale_y= 40;
 %scale barx
 hold on;x1= (1250*srF)-(100*srF);x2=1250*srF;p1=plot([x1 x2],[ov_min-10 ov_min-10],'-','Color','k','LineWidth',1.5);
 %scale bary
 hold on;y2= (ov_min-10)+scale_y;y1=(ov_min-10);p2=plot([x2 x2],[y1 y2],'-','Color','k','LineWidth',1.5); 
title([num2str(data(temp(i)).patching_date) data(temp(i)).celID ]);
else
plot(1,1);set(gca,'box','off');axis off;
end
end

elseif ty==2 %Rheobase
 trace_1xRheo_2ap=[];
for  i=1:length(find(idx==1))
    temp=find(idx==1);
    if isempty(data(temp(i)).Rheobase)==0;
        %find Rheobase 
        rheo(i)=data(temp(i)).Rheobase.rheo;
        idx_1xRheo(i)=find(data(temp(i)).Rheobase.stimvec==1*data(temp(i)).Rheobase.rheo);
        trace_1xRheo(:,i)=data(temp(i)).Rheobase.traces(:,idx_1xRheo(i));
         if isempty(find(data(temp(i)).Rheobase.spikecount>=2))==0
    md=find(data(temp(i)).Rheobase.spikecount>=2);
    trace_2ap(:,i)=data(temp(i)).Rheobase.traces(:,md(1));
   else
      md=find(data(temp(i)).Rheobase.spikecount>=1);
    trace_2ap(:,i)=data(temp(i)).Rheobase.traces(:,md(1));
   end
        trace_1xRheo_2ap(:,i)=data(temp(i)).Rheobase.traces(:,idx_1xRheo(i)+1);
        sp_time_1x=spike_times(trace_1xRheo(:,i),0.1);
        if idx_1xRheo(i)<9
       % idx_2xRheo(i)=find(data(temp(i)).Rheobase.stimvec==data(temp(i)).Rheobase.rheo(end));
       idx_2xRheo(i)=16;
        else
        %idx_2xRheo(i)=find(data(temp(i)).Rheobase.stimvec==data(temp(i)).Rheobase.stimvec(end));  
        idx_2xRheo(i)=16;
        end
        trace_2xRheo(:,i)=data(temp(i)).Rheobase.traces(:,idx_2xRheo(i));
    else
        rheo(i)=NaN;
        idx_1xRheo(i)=NaN;
        trace_2ap(:,i)=ones(length(data(temp(1)).Rheobase.traces),1)*NaN;
        trace_2xRheo(:,i)=ones(length(data(temp(1)).Rheobase.traces),1)*NaN;
    end
    sp_time_2x=spike_times(trace_2xRheo(:,i),0.1);
    try
    F2xRheo(:,i)=1/((sp_time_2x(2)-sp_time_2x(1))/sr);
    catch
    F2xRheo(:,i)=0;
    end
end 
fig1=figure;set(fig1, 'Position', [200, -200, 800, 800]);set(gcf,'color','w');
for i=1:length(find(idx==1))
    subplot(siz,siz,i);
   plot(trace_2ap(:,i),'Color',col,'LineWidth',1); box off;
  hold on; plot(trace_2xRheo(:,i),'Color','r','LineWidth',1); box off;
end

else
    
end

end