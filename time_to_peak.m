function [time_tp]=time_to_peak(trac,bs,act,sr,cc)

st_bs=3.5*std(diff(trac(bs)));
if cc==0
temp=find(diff(trac(act))<-st_bs);
figure;plot(diff(trac(act)));hold on;plot([0 100],[-st_bs -st_bs],'--k');
if isempty(temp)==0
time_tp=temp(1);
else
    time_tp=NaN;
end

else cc==1
temp=find(diff(trac(act))>st_bs);
figure;plot(diff(trac(act)));hold on;plot([0 100],[st_bs st_bs],'--k');
if isempty(temp)==0
time_tp=temp(1);
else
    time_tp=NaN;
end
end
end