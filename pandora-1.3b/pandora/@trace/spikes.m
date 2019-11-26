function obj = spikes(t, a_period, plotit, minamp)

% spikes - Convert trace to spikes object for spike timing calculations.
%
% Usage:
% obj = spikes(trace, a_period, plotit, minamp)
%
%   Parameters:
%	trace: A trace object.
%	a_period: A period object denoting the part of trace of interest 
%		(optional, if empty vector, taken as wholePeriod).
%	plotit: If non-zero, a plot is generated for showing spikes found
%		(optional).
%	minamp: minimum amplitude that must be reached if using findFilteredSpikes.
%		--> adjust as needed to discriminate spikes from EPSPs.
%		(optional)
%
% Description:
%   Creates a spikes object.
%		
% See also: spikes
%
% $Id: spikes.m 1335 2012-04-19 18:04:32Z cengique $
%
% Author: 
%   Cengiz Gunay <cgunay@emory.edu>, 2004/07/30
% Modified:
% - added minamp parameter. J. Edgerton, 10/10/2005

% Copyright (c) 2007 Cengiz Gunay <cengique@users.sf.net>.
% This work is licensed under the Academic Free License ("AFL")
% v. 3.0. To view a copy of this license, please look at the COPYING
% file distributed with this software or visit
% http://opensource.org/licenses/afl-3.0.php.

if nargin == 0 % Called with no params
  error('Need trace parameter.');
end

if ~ exist('a_period', 'var') || isempty(a_period)
  a_period = periodWhole(t);
end
  
if ~ exist('plotit', 'var')
  plotit = 0;
end

if ~ exist('minamp', 'var')
  if isfield(t.props, 'threshold')
    minamp = t.props.threshold;
  else
    minamp = [];    % default is in findFilteredSpikes.
  end
end

% if the trace object contains multiple traces, return array of spikes
% objects 
num_traces = size(t.data, 2);
if num_traces > 1
  obj = repmat(spikes, 1, num_traces);
  for trace_num = 1:num_traces
    obj(trace_num) = ...
        spikes(set(t, 'data', t.data(:, trace_num)), a_period, plotit, minamp);
  end
  return
end

if (a_period.end_time - a_period.start_time) > 0

% Choose an appropriate spike finder here and indicate in id.
if isfield(t.props, 'spike_finder') && ...
      t.props.spike_finder == 2 && ...
      isfield(t.props, 'threshold')
  % Scale to mV for spike finder
  mV_factor = 1e3 * t.dy;
  
  if plotit > 0, plot_str = {'plot'}; else plot_str = {}; end

  % Li Su's findspikes requires a hack to pass it 1kHz sampling rate to
  % return spike times in terms of dt
  [times, peaks, n] = ...
      findspikes(t.data(a_period.start_time:a_period.end_time) * mV_factor, ...
		 1, t.props.threshold, plot_str{:});
elseif isfield(t.props, 'spike_finder') && ...
      t.props.spike_finder == 3 && ...
      isfield(t.props, 'threshold') 
    % Scale to mV for spike finder
    mV_factor = 1e3 * t.dy;
    
    % Alfonso's old findspikes. Simple, but works
    [times, peaks, n] = ...
      findspikes_old(t.data(a_period.start_time:a_period.end_time) * mV_factor, ...
		 t.props.threshold, plotit);
else 
  % Assume spike_finder == 1 for filtered method.
  % Pass t.props:
  [times, peaks, n] = ...
      findFilteredSpikes(t, a_period, plotit, minamp, t.props);
end

obj = spikes(times, length(t.data), t.dt, t.id);

else

obj = spikes([], length(t.data), t.dt, t.id);

end
