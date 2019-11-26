function h = plot(t, title_str)

% plot - Plots spikes.
%
% Usage: 
% h = plot(t)
%
% Description:
%
%   Parameters:
%	t: A spikes object.
%	title_str: (Optional) String to append to plot title.
%
%   Returns:
%	h: Handle to figure object.
%
% See also: spikes, plot_abstract
%
% $Id: plot.m 1335 2012-04-19 18:04:32Z cengique $
%
% Author: Cengiz Gunay <cgunay@emory.edu>, 2004/08/04

% Copyright (c) 2007 Cengiz Gunay <cengique@users.sf.net>.
% This work is licensed under the Academic Free License ("AFL")
% v. 3.0. To view a copy of this license, please look at the COPYING
% file distributed with this software or visit
% http://opensource.org/licenses/afl-3.0.php.

if ~ exist('title_str', 'var')
  title_str = '';
end

s = size(t);
if max(s) > 1
  % Column vector
  if s(1) > s(2)
    % Make a vertical stack plot (default)
    orientation = 'y';
  else
    orientation = 'x';		% or horizontal
  end
  plotFigure(plot_stack(num2cell(plotData(t)), [], orientation, title_str));
else
  h = plotFigure(plotData(t, title_str));
end