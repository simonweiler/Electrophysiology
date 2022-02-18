function [time_tp trace_smooth]=time_to_peak(trac,bs,act,sr,cc,fc)

%st_bs=4*std(diff(trac(bs)));
filter_used='gaussian';
bs_diff=smoothdata(diff(trac(bs)),filter_used);
act_diff=smoothdata(diff(trac(act)),filter_used);
st_bs=fc*std(bs_diff);
trace_smooth=trac(act);

if cc==0
%temp=find(diff(trac(act))<-st_bs);
temp=find(act_diff<-st_bs);
figure;plot(act_diff);hold on;plot([0 1000],[-st_bs -st_bs],'--k');
if isempty(temp)==0
time_tp=temp(1)/sr;
else
    time_tp=NaN;
end

else cc==1
%temp=find(diff(trac(act))>st_bs);
temp=find(act_diff>st_bs);
figure;plot(act_diff);hold on;plot([0 1000],[-st_bs -st_bs],'--k');
if isempty(temp)==0
time_tp=temp(1)/sr;
else
    time_tp=NaN;
end

end
end