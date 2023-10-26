function [time_tp trace_smooth]=time_to_peak_callosal(trac,bs,act,sr,cc,fc,cellnr)


filter_used='gaussian';
filter_trac=smoothdata(trac,filter_used,100);
bs_diff=diff(filter_trac(bs));
act_diff=diff(filter_trac(act));
st_bs=fc*std(bs_diff);
trace_smooth=trac(act);


if cc==0
%temp=find(diff(trac(act))<-st_bs);
temp=find(act_diff<-st_bs);

if isempty(temp)==0 & max(abs(filter_trac(act))>30)==1
     %if temp(1)>20
time_tp=temp(1)/sr;
     %else
   %idx=find(temp>=max(diff(temp)));
  %time_tp=temp(idx(1))/sr;
   %  end
else
    time_tp=NaN;
end
fig7= figure;set(fig7, 'Position', [400, 500, 400, 400]);set(gcf,'color','w');
subplot(2,1,1);plot(filter_trac(act),'Color','k');title([num2str(cellnr) ' ' num2str(time_tp)]);box off;
subplot(2,1,2);plot(act_diff,'Color','r');hold on;plot([0 1000],[-st_bs -st_bs],'--k');box off;

else cc==1
%temp=find(diff(trac(act))>st_bs);
temp=find(act_diff>st_bs);
if isempty(temp)==0 & max(abs(filter_trac(act))>30)==1
    %if temp(1)>20
time_tp=temp(1)/sr;
   % else
  % idx=find(temp>=max(diff(temp)));
  % time_tp=temp(idx(1))/sr;
    %end
else
    time_tp=NaN;
end

fig7= figure;set(fig7, 'Position', [400, 500, 400, 400]);set(gcf,'color','w');
subplot(2,1,1);plot(filter_trac(act),'Color','k');title([num2str(cellnr) ' ' num2str(time_tp)]);box off;
subplot(2,1,2);plot(act_diff,'Color','b');hold on;plot([0 1000],[st_bs st_bs],'--k');box off;


end

end