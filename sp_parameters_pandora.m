function  [active_p] = sp_parameters_pandora(input,nr)

%this will use PANDORA TOOLBOX https://github.com/cengique/pandora-matlab
%input=trace_pyr is the spike trace (trace at Rheobase is recommended)
%nr= if trace has more than one spike use second spike in the trace 
%
%Parameters are:
% %%%%[ 1]    'MinVm'        
%     [ 2]    'PeakVm'       
%     [ 3]    'InitVm'       
%     [ 4]    'InitVmBySlope'
%     [ 5]    'MaxVmSlope'   
%     [ 6]    'HalfVm'       
%     [ 7]    'Amplitude'    
%     [ 8]    'MaxAHP'       
%     [ 9]    'DAHPMag'      
%     [10]    'InitTime'     
%     [11]    'RiseTime'     
%     [12]    'FallTime'     
%     [13]    'MinTime'      
%     [14]    'BaseWidth'    
%     [15]    'HalfWidth'    
%     [16]    'FixVWidth'    
%     [17]    'Index'        
%     [18]    'Time'
%Pandora parameters
% https://github.com/cengique/pandora-matlab/blob/master/doc/tutorials/incf/extracting-spike-info.markdown
% 
% The basic description of the main measures are:
% 
% `SpikeAmplitude`	Amplitude of spike from threshold point to tip in mV.
% `SpikeBaseWidth`	Width in ms at the threshold level of the spike.
% `SpikeDAHPMag`	(Experimental) Double afterhyperpolarization (AHP) magnitude in mV.
% `SpikeFallTime`	Time of repolarization in ms from the tip of the spike to threshold crossing.
% `SpikeFixVWidth`	Width in ms of the spike at a fixed voltage value (default: -10 mV, modify in `spike_shape/getResults.m`)
% `SpikeHalfVm`	Voltage in mV at spike half-height.
% `SpikeHalfWidth`	Width in ms at spike half-height.
% `SpikeInitTime`	Spike initiation time (at threshold) in ms from start of extracted spike shape.
% `SpikeInitVmBySlope`	Spike initiation point (threshold) voltage in mV estimated by slope crossing method (see `spike_shape/calcInitVm.m` for different methods)
% `SpikeInitVm`	Spike threshold voltage in mV estimated by *best* method (First maximum curvature, then slope method is tried. See `spike_shape/calcInitVm.m`).
% `SpikeMaxAHP`	Maximal magnitude in mV of afterhyperpolarization (AHP).
% `SpikeMaxVmSlope`	Maximal slope of voltage in mV/ms during spike depolarization (rise).
% `SpikeMinTime`	Minimum of spike shape in ms starting from beginning of extracted spike shape.
% `SpikeMinVm`	Absolute voltage in mV at spike minimum.
% `SpikePeakVm`	Absolute voltage in mV at spike maximum (peak).
% `SpikeRiseTime`	Time of depolarization in ms from threshold point to the spike tip.
% `Spikes`	Number of spikes found.

for i=1:size(input,2);
trace_curr=input(:,i);
dt = 1e-4;dy = 1e-3;
props = struct('spike_finder', 2, 'threshold', 0);
traces_analysis=trace(trace_curr,dt, dy, 'Analysis', props);
alltrace_info = getProfileAllSpikes(traces_analysis);
parameters=[];
parameters=alltrace_info.spikes_db.data; 
try
active_p(:,i)=parameters(nr,:);
catch
active_p(:,i)=parameters(1,:);
end;
end
end