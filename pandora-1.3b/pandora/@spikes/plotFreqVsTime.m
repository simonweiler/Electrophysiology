function a_plot = plotFreqVsTime(s, title_str, props)

% plotFreqVsTime - Plots a frequency-time graph from the spikes object.
%
% Usage: 
% a_plot = plotFreqVsTime(s, title_str, props)
%
% Description:
%   If s is a vector of spikes objects, returns a vector of plot objects.
%
%   Parameters:
%	s: A spikes object.
%	title_str: (Optional) String to append to plot title.
%	props: A structure with any optional properties.
%	  timeScale: 's' for seconds, or 'ms' for milliseconds.
%	  type: If 'simple' plots 1/is for each spike time, 
%		'manhattan' uses flat lines of 1/isi height between spike times (default).
%	  (others passed to plot_abstract)
%
%   Returns:
%	a_plot: A plot_abstract object that can be visualized.
%	title_str: (Optional) String to append to plot title.
%
% See also: trace, plot_abstract
%
% $Id: plotFreqVsTime.m 1335 2012-04-19 18:04:32Z cengique $
%
% Author: Cengiz Gunay <cgunay@emory.edu>, 2006/05/05

% Copyright (c) 2007 Cengiz Gunay <cengique@users.sf.net>.
% This work is licensed under the Academic Free License ("AFL")
% v. 3.0. To view a copy of this license, please look at the COPYING
% file distributed with this software or visit
% http://opensource.org/licenses/afl-3.0.php.

if ~ exist('props', 'var')
  props = struct;
end

if ~ exist('title_str', 'var')
  title_str = '';
end

% If input is an array, then return array of plots
num_dbs = length(s);
if num_dbs > 1 
  % Create array of plots
  [a_plot(1:num_dbs)] = deal(plot_abstract);
  for plot_num = 1:num_dbs
    a_plot(plot_num) = plotFreqVsTime(s(plot_num), title_str);
  end
  return;
end

if ~isfield(props, 'timeScale')
  props.timeScale = 'ms';
end

switch props.timeScale
  case 's'
    time_factor = s.dt;
    x_label = 'time [s]';
  case 'ms'
    time_factor = s.dt * 1e3;
    x_label = 'time [ms]';
end

% Remove all '_' characters, because they interfere with TeX interpretation
class_name = strrep(class(s), '_', ' ');

freqs = 1 ./ getISIs(s) ./ s.dt;

if ~ isfield(props, 'type')
  props.type = 'manhattan';
end

switch props.type
  case 'simple'
    if length(s.times) > 1
      stimes = [(s.times(1) - .1 * time_factor) s.times(1:end-1) ...
		(s.times(end - 1) + .1 * time_factor) ] * time_factor;
    else
      stimes = [ 0 0 ];
    end
    freqs = [ 0 freqs 0 ];

  case 'manhattan'
    if length(s.times) > 1
      stimes = s.times * time_factor;
      num_stimes = 2 * length(stimes) + 1;
      new_stimes = zeros(1, num_stimes);

      % duplicate values for drawing ISIs as flat lines
      new_stimes(2:2:num_stimes) = stimes;
      new_stimes(3:2:num_stimes) = stimes;
      stimes = new_stimes;

      % freqs
      new_freqs = zeros(1, num_stimes);
      new_freqs(3:2:(num_stimes-1)) = freqs;
      new_freqs(4:2:num_stimes) = freqs;
      freqs = new_freqs;
      
    else
      stimes = [ 0 0 ];
      freqs = [ 0 0 ];
    end
  
  otherwise
    error(['Error: Plot type ' props.type ' not known.']);
end

class_name = strrep(class(s), '_', ' ');
the_title = [ sprintf('%s freq-vs-time: %s', class_name, s.id) title_str ...
            ];
the_legend = s.id;
if isfield(props, 'quiet')
  if ~ isempty(title_str)
    the_title = title_str;
    the_legend = title_str;
  end
else
end

a_plot = plot_abstract({stimes, freqs}, ...
		       {x_label, 'firing rate [Hz]'}, ...
		       properTeXLabel(the_title), { properTeXLabel(the_legend) }, 'plot', props);
