%Plot all cells to assess the quality of recording 
function plot_slice_opto(Ephys, cell_idx, sol, sample_rate)
%% OUPTPUT
    %plots PSCs/PSPs of different freq of light stimulations for a given cell 

%% INPUTS
    %Ephys: slice opto data structure
    %cell_idx=cells ID of which PSCs/PSPs to plot
    %sol: Solution used (K-gluc=1;Cs-gluc=2)
    %sample_rate: rate in Hz used for sampling data (e.g. 20000)

%% Set up parameters 
plot_time=[1/sample_rate:1/sample_rate:6];
laser_on=[0.25 1.25 2.25 3.25 4.25];
laser_off=[0.35 1.35 2.35 3.35 4.35];
laser_on2=[0.25:0.2:5.2];
laser_off2=[0.26:0.2:5.21];
laser_on3=[0.25:0.1:5.2];
laser_off3=[0.26:0.1:5.21];
%% PLOT STARTS HERE

fig7=figure;set(gcf,'color','w');set(fig7, 'Position', [0, 200 ,1200, 700]);

%K-gluconate so therefore both EPSC and/or EPSP
if sol==1
     %% 1HZ so 5 pulses, 100 ms
     %plot EPSC or EPSP 1Hz, 5 pulses 100 ms
    if ~isempty(Ephys(cell_idx).sub_traces_train)==1
      for i=1:length(Ephys(cell_idx).clamp1)
        %EPSC plot
        if Ephys(cell_idx).clamp1(i)==0
        try
        subplot(3,2,1);p1=plot(plot_time,Ephys(cell_idx).sub_traces_train(1:length(plot_time),i),'Color','k','Linewidth',0.5);p1.Color(4) = 0.8;
        catch
        subplot(3,2,1);p1=plot(plot_time(1:length(Ephys(cell_idx).sub_traces_train)),Ephys(cell_idx).sub_traces_train(:,i),'Color','k','Linewidth',0.5);p1.Color(4) = 0.8;
        end
        hold on;line(get(gca, 'xlim'), [0 0], 'color', 'r', 'linestyle', ':','LineWidth',1.5);box off;hold on;
            %plot laser on and off as lines
            for k=1:length(laser_on)
            p1=plot([laser_on(k), laser_on(k)], [get(gca, 'ylim')], 'color', 'b', 'linestyle', ':','LineWidth',0.8);p1.Color(4) = 0.8;hold on;
            p1=plot([laser_off(k), laser_off(k)], [get(gca, 'ylim')], 'color', 'b', 'linestyle', ':','LineWidth',0.8);p1.Color(4) = 0.8;hold on;
            end
        temp=get(gca, 'ylim');ylim([temp(1)-50 temp(2)]);
        hold on;line([5.8, 6], [temp(1)-50 temp(1)-50], 'color', 'k', 'linestyle', '-','LineWidth',1);hold on;
        hold on;line([5.8, 5.8], [temp(1)-50 temp(1)-50+40], 'color', 'k', 'linestyle', '-','LineWidth',1);hold on;axis off;hold on;
        title([Ephys(cell_idx).animal_name '  ' num2str(Ephys(cell_idx).patching_date) '  ' Ephys(cell_idx).celID ' cell '  num2str(cell_idx) ' opto ' num2str(Ephys(cell_idx).optovariant)...
            ' pair ' num2str(Ephys(cell_idx).pair) ' EPSC']);
        %EPSP plot
        else
        subplot(3,2,2);p1=plot(plot_time,Ephys(cell_idx).sub_traces_train(1:length(plot_time),i),'Color','k','Linewidth',0.5);p1.Color(4) = 0.8;
        hold on;line(get(gca, 'xlim'), [0 0], 'color', 'r', 'linestyle', ':','LineWidth',1.5);box off;
            %plot laser on and off as lines
            for k=1:length(laser_on)
            p1=plot([laser_on(k), laser_on(k)], [get(gca, 'ylim')], 'color', 'b', 'linestyle', ':','LineWidth',0.8);p1.Color(4) = 0.8;hold on;
            p1=plot([laser_off(k), laser_off(k)], [get(gca, 'ylim')], 'color', 'b', 'linestyle', ':','LineWidth',0.8);p1.Color(4) = 0.8;hold on;
            end
        temp=get(gca, 'ylim');ylim([temp(1) temp(2)+5]);
        hold on;line([5.8, 6], [temp(1) temp(1)], 'color', 'k', 'linestyle', '-','LineWidth',1);hold on;
        hold on;line([5.8, 5.8], [temp(1) temp(1)+2], 'color', 'k', 'linestyle', '-','LineWidth',1);hold on;axis off;
        title('EPSP');
        end
     end
    end
%% 5 Hz recording so 25 pulses, 10 ms pulses
   if ~isempty(Ephys(cell_idx).sub_traces_high)==1
     for i=1:length(Ephys(cell_idx).clamp3)
        if Ephys(cell_idx).clamp3(i)==0
        subplot(3,2,3);p1=plot(plot_time,Ephys(cell_idx).sub_traces_high(1:length(plot_time),i),'Color','k','Linewidth',0.5);p1.Color(4) = 0.8;
        hold on;line(get(gca, 'xlim'), [0 0], 'color', 'r', 'linestyle', ':','LineWidth',1.5);box off;hold on;
            %plot laser on and off as lines
            for k=1:length(laser_on2)
            p1=plot([laser_on2(k), laser_on2(k)], [get(gca, 'ylim')], 'color', 'b', 'linestyle', ':','LineWidth',0.8);p1.Color(4) = 0.8;hold on;
            %p1=plot([laser_off2(k), laser_off2(k)], [get(gca, 'ylim')], 'color', 'b', 'linestyle', ':','LineWidth',0.8);p1.Color(4) = 0.8;hold on;
            end
        temp=get(gca, 'ylim');ylim([temp(1)-50 temp(2)]);
        hold on;line([5.8, 6], [temp(1)-50 temp(1)-50], 'color', 'k', 'linestyle', '-','LineWidth',1);hold on;
        hold on;line([5.8, 5.8], [temp(1)-50 temp(1)-50+40], 'color', 'k', 'linestyle', '-','LineWidth',1);hold on;axis off;hold on;
        else
        subplot(3,2,4);p1=plot(plot_time,Ephys(cell_idx).sub_traces_high(1:length(plot_time),i),'Color','k','Linewidth',0.5);p1.Color(4) = 0.8;
        hold on;line(get(gca, 'xlim'), [0 0], 'color', 'r', 'linestyle', ':','LineWidth',1.5);box off;
            %plot laser on and off as lines
            for k=1:length(laser_on2)
            p1=plot([laser_on2(k), laser_on2(k)], [get(gca, 'ylim')], 'color', 'b', 'linestyle', ':','LineWidth',0.8);p1.Color(4) = 0.8;hold on;
            %p1=plot([laser_off2(k), laser_off2(k)], [get(gca, 'ylim')], 'color', 'b', 'linestyle', ':','LineWidth',0.8);p1.Color(4) = 0.8;hold on;
            end
         temp=get(gca, 'ylim');ylim([temp(1) temp(2)+5]);
         hold on;line([5.8, 6], [temp(1) temp(1)], 'color', 'k', 'linestyle', '-','LineWidth',1);hold on;
         hold on;line([5.8, 5.8], [temp(1) temp(1)+2], 'color', 'k', 'linestyle', '-','LineWidth',1);hold on;axis off;
        end
     end
    end
%% %% 10 Hz recording so 50 pulses, 10 ms pulses
   if ~isempty(Ephys(cell_idx).sub_traces_high)==1
     for i=1:length(Ephys(cell_idx).clamp4)
        if Ephys(cell_idx).clamp4(i)==0
        subplot(3,2,5);p1=plot(plot_time,Ephys(cell_idx).sub_traces_highf(1:length(plot_time),i),'Color','k','Linewidth',0.5);p1.Color(4) = 0.8;
        hold on;line(get(gca, 'xlim'), [0 0], 'color', 'r', 'linestyle', ':','LineWidth',1.5);box off;hold on;
            %plot laser on and off as lines
            for k=1:length(laser_on3)
            p1=plot([laser_on3(k), laser_on3(k)], [get(gca, 'ylim')], 'color', 'b', 'linestyle', ':','LineWidth',0.8);p1.Color(4) = 0.8;hold on;
            %p1=plot([laser_off3(k), laser_off3(k)], [get(gca, 'ylim')], 'color', 'b', 'linestyle', ':','LineWidth',0.8);p1.Color(4) = 0.8;hold on;
            end
        temp=get(gca, 'ylim');ylim([temp(1)-50 temp(2)]);
        hold on;line([5.8, 6], [temp(1)-50 temp(1)-50], 'color', 'k', 'linestyle', '-','LineWidth',1);hold on;
        hold on;line([5.8, 5.8], [temp(1)-50 temp(1)-50+40], 'color', 'k', 'linestyle', '-','LineWidth',1);hold on;axis off;hold on;
        else
        subplot(3,2,6);p1=plot(plot_time,Ephys(cell_idx).sub_traces_highf(1:length(plot_time),i),'Color','k','Linewidth',0.5);p1.Color(4) = 0.8;
        hold on;line(get(gca, 'xlim'), [0 0], 'color', 'r', 'linestyle', ':','LineWidth',1.5);box off;
            %plot laser on and off as lines
            for k=1:length(laser_on3)
            p1=plot([laser_on3(k), laser_on3(k)], [get(gca, 'ylim')], 'color', 'b', 'linestyle', ':','LineWidth',0.8);p1.Color(4) = 0.8;hold on;
            %p1=plot([laser_off3(k), laser_off3(k)], [get(gca, 'ylim')], 'color', 'b', 'linestyle', ':','LineWidth',0.8);p1.Color(4) = 0.8;hold on;
            end
         temp=get(gca, 'ylim');ylim([temp(1) temp(2)+5]);
         hold on;line([5.8, 6], [temp(1) temp(1)], 'color', 'k', 'linestyle', '-','LineWidth',1);hold on;
         hold on;line([5.8, 5.8], [temp(1) temp(1)+2], 'color', 'k', 'linestyle', '-','LineWidth',1);hold on;axis off;
        end
     end
    end

%CS-gluc SOLUTION so EPSC and IPSC in the same cell 
else sol==2
    %%      %% 1HZ so 5 pulses, 100 ms
     %plot EPSC or EPSP 1Hz, 5 pulses 100 ms
    if ~isempty(Ephys(cell_idx).sub_traces_train)==1
      for i=1:length(Ephys(cell_idx).ipsc_log)
        %EPSC plot
        if Ephys(cell_idx).ipsc_log(i)==0
        subplot(3,2,1);p1=plot(plot_time,Ephys(cell_idx).sub_traces_train(1:length(plot_time),i),'Color','k','Linewidth',0.5);p1.Color(4) = 0.8;
        hold on;line(get(gca, 'xlim'), [0 0], 'color', 'r', 'linestyle', ':','LineWidth',1.5);box off;hold on;
            %plot laser on and off as lines
            for k=1:length(laser_on)
            p1=plot([laser_on(k), laser_on(k)], [get(gca, 'ylim')], 'color', 'b', 'linestyle', ':','LineWidth',0.8);p1.Color(4) = 0.8;hold on;
            p1=plot([laser_off(k), laser_off(k)], [get(gca, 'ylim')], 'color', 'b', 'linestyle', ':','LineWidth',0.8);p1.Color(4) = 0.8;hold on;
            end
        temp=get(gca, 'ylim');ylim([temp(1)-50 temp(2)]);
        hold on;line([5.8, 6], [temp(1)-50 temp(1)-50], 'color', 'k', 'linestyle', '-','LineWidth',1);hold on;
        hold on;line([5.8, 5.8], [temp(1)-50 temp(1)-50+40], 'color', 'k', 'linestyle', '-','LineWidth',1);hold on;axis off;hold on;
         title([Ephys(cell_idx).animal_name '  ' num2str(Ephys(cell_idx).patching_date) '  ' Ephys(cell_idx).celID ' cell '  num2str(cell_idx) ' opto ' num2str(Ephys(cell_idx).optovariant)...
            ' pair ' num2str(Ephys(cell_idx).pair) ' EPSC']);
        %IPSC plot
        else
        subplot(3,2,2);p1=plot(plot_time,Ephys(cell_idx).sub_traces_train(1:length(plot_time),i),'Color','k','Linewidth',0.5);p1.Color(4) = 0.8;
        hold on;line(get(gca, 'xlim'), [0 0], 'color', 'r', 'linestyle', ':','LineWidth',1.5);box off;
            %plot laser on and off as lines
            for k=1:length(laser_on)
            p1=plot([laser_on(k), laser_on(k)], [get(gca, 'ylim')], 'color', 'b', 'linestyle', ':','LineWidth',0.8);p1.Color(4) = 0.8;hold on;
            p1=plot([laser_off(k), laser_off(k)], [get(gca, 'ylim')], 'color', 'b', 'linestyle', ':','LineWidth',0.8);p1.Color(4) = 0.8;hold on;
            end
        temp=get(gca, 'ylim');ylim([temp(1) temp(2)+5]);
        hold on;line([5.8, 6], [temp(1) temp(1)], 'color', 'k', 'linestyle', '-','LineWidth',1);hold on;
        hold on;line([5.8, 5.8], [temp(1) temp(1)+50], 'color', 'k', 'linestyle', '-','LineWidth',1);hold on;axis off;
        title('IPSC');
        end
     end
    end
%%   plot EPSC or IPSC 5Hz, 25 pulses 10 ms
    if ~isempty(Ephys(cell_idx).sub_traces_train)==1
      for i=1:length(Ephys(cell_idx).ipsc_log_high)
        %EPSC plot
        if Ephys(cell_idx).ipsc_log_high(i)==0
        subplot(3,2,3);p1=plot(plot_time,Ephys(cell_idx).sub_traces_high(1:length(plot_time),i),'Color','k','Linewidth',0.5);p1.Color(4) = 0.8;
        hold on;line(get(gca, 'xlim'), [0 0], 'color', 'r', 'linestyle', ':','LineWidth',1.5);box off;hold on;
            %plot laser on and off as lines
            for k=1:length(laser_on2)
            p1=plot([laser_on2(k), laser_on2(k)], [get(gca, 'ylim')], 'color', 'b', 'linestyle', ':','LineWidth',0.8);p1.Color(4) = 0.8;hold on;
            %p1=plot([laser_off2(k), laser_off2(k)], [get(gca, 'ylim')], 'color', 'b', 'linestyle', ':','LineWidth',0.8);p1.Color(4) = 0.8;hold on;
            end
        temp=get(gca, 'ylim');ylim([temp(1)-50 temp(2)]);
        hold on;line([5.8, 6], [temp(1)-50 temp(1)-50], 'color', 'k', 'linestyle', '-','LineWidth',1);hold on;
        hold on;line([5.8, 5.8], [temp(1)-50 temp(1)-50+40], 'color', 'k', 'linestyle', '-','LineWidth',1);hold on;axis off;hold on;
        title([Ephys(cell_idx).animal_name '  ' num2str(Ephys(cell_idx).patching_date) '  ' Ephys(cell_idx).celID ' cell '  num2str(cell_idx) ' EPSC']);
        %IPSC plot
        else
        subplot(3,2,4);p1=plot(plot_time,Ephys(cell_idx).sub_traces_high(1:length(plot_time),i),'Color','k','Linewidth',0.5);p1.Color(4) = 0.8;
        hold on;line(get(gca, 'xlim'), [0 0], 'color', 'r', 'linestyle', ':','LineWidth',1.5);box off;
            %plot laser on and off as lines
            for k=1:length(laser_on2)
            p1=plot([laser_on2(k), laser_on2(k)], [get(gca, 'ylim')], 'color', 'b', 'linestyle', ':','LineWidth',0.8);p1.Color(4) = 0.8;hold on;
            %p1=plot([laser_off2(k), laser_off2(k)], [get(gca, 'ylim')], 'color', 'b', 'linestyle', ':','LineWidth',0.8);p1.Color(4) = 0.8;hold on;
            end
        temp=get(gca, 'ylim');ylim([temp(1) temp(2)+5]);
        hold on;line([5.8, 6], [temp(1) temp(1)], 'color', 'k', 'linestyle', '-','LineWidth',1);hold on;
        hold on;line([5.8, 5.8], [temp(1) temp(1)+50], 'color', 'k', 'linestyle', '-','LineWidth',1);hold on;axis off;
        title('IPSC');
        end
     end
    end
%% plot EPSC or IPSC 10Hz, 50 pulses 10 ms
    if ~isempty(Ephys(cell_idx).sub_traces_train)==1
      for i=1:length(Ephys(cell_idx).ipsc_log_highf)
        %EPSC plot
        if Ephys(cell_idx).ipsc_log_highf(i)==0
        subplot(3,2,5);p1=plot(plot_time,Ephys(cell_idx).sub_traces_highf(1:length(plot_time),i),'Color','k','Linewidth',0.5);p1.Color(4) = 0.8;
        hold on;line(get(gca, 'xlim'), [0 0], 'color', 'r', 'linestyle', ':','LineWidth',1.5);box off;hold on;
            %plot laser on and off as lines
            for k=1:length(laser_on3)
            p1=plot([laser_on3(k), laser_on3(k)], [get(gca, 'ylim')], 'color', 'b', 'linestyle', ':','LineWidth',0.8);p1.Color(4) = 0.8;hold on;
            %p1=plot([laser_off3(k), laser_off3(k)], [get(gca, 'ylim')], 'color', 'b', 'linestyle', ':','LineWidth',0.8);p1.Color(4) = 0.8;hold on;
            end
        temp=get(gca, 'ylim');ylim([temp(1)-50 temp(2)]);
        hold on;line([5.8, 6], [temp(1)-50 temp(1)-50], 'color', 'k', 'linestyle', '-','LineWidth',1);hold on;
        hold on;line([5.8, 5.8], [temp(1)-50 temp(1)-50+40], 'color', 'k', 'linestyle', '-','LineWidth',1);hold on;axis off;hold on;
        title([Ephys(cell_idx).animal_name '  ' num2str(Ephys(cell_idx).patching_date) '  ' Ephys(cell_idx).celID ' cell '  num2str(cell_idx) ' EPSC']);
        %IPSC plot
        else
        subplot(3,2,6);p1=plot(plot_time,Ephys(cell_idx).sub_traces_highf(1:length(plot_time),i),'Color','k','Linewidth',0.5);p1.Color(4) = 0.8;
        hold on;line(get(gca, 'xlim'), [0 0], 'color', 'r', 'linestyle', ':','LineWidth',1.5);box off;
            %plot laser on and off as lines
            for k=1:length(laser_on3)
            p1=plot([laser_on3(k), laser_on3(k)], [get(gca, 'ylim')], 'color', 'b', 'linestyle', ':','LineWidth',0.8);p1.Color(4) = 0.8;hold on;
            %p1=plot([laser_off3(k), laser_off3(k)], [get(gca, 'ylim')], 'color', 'b', 'linestyle', ':','LineWidth',0.8);p1.Color(4) = 0.8;hold on;
            end
        temp=get(gca, 'ylim');ylim([temp(1) temp(2)+5]);
        hold on;line([5.8, 6], [temp(1) temp(1)], 'color', 'k', 'linestyle', '-','LineWidth',1);hold on;
        hold on;line([5.8, 5.8], [temp(1) temp(1)+50], 'color', 'k', 'linestyle', '-','LineWidth',1);hold on;axis off;
        title('IPSC');
        end
     end
    end
end
end
