function a_plot = plotrows(a_tests_db, axis_limits, orientation, title_str, props)

% plotrows - Creates a plot_stack describing the db rows.
%
% Usage:
% a_plot = plotrows(a_tests_db, axis_limits, orientation, props)
%
% Description:
%
%   Parameters:
%	a_tests_db: A tests_db object.
%	axis_limits: If given, all plots contained will have these axis limits.
%	orientation: Stack orientation 'x' for horizontal, 'y' for vertical, etc.
%	title_str: Optional title string.
%	props: A structure with any optional properties passed to plot_stack.
%		
%   Returns:
%	a_plot: A plot_stack object that can be plotted.
%
% See also: plot_abstract, plotFigure
%
% $Id: plotrows.m 1335 2012-04-19 18:04:32Z cengique $
%
% Author: Cengiz Gunay <cgunay@emory.edu>, 2004/11/09

% Copyright (c) 2007 Cengiz Gunay <cengique@users.sf.net>.
% This work is licensed under the Academic Free License ("AFL")
% v. 3.0. To view a copy of this license, please look at the COPYING
% file distributed with this software or visit
% http://opensource.org/licenses/afl-3.0.php.

if ~ exist('props', 'var')
  props = struct([]);
end

if ~ exist('axis_limits', 'var')
  axis_limits = [];
end

if ~ exist('orientation', 'var')
  orientation = 'y';
end

if ~ exist('title_str', 'var')
  title_str = '';
end

num_rows = dbsize(a_tests_db, 1);
plots = cell(num_rows, 1);
%border = [0 0 0 0];
for row_num=1:num_rows
  % inverse order for plot_stack
  if row_num == 1
    row_plot = ...
	plotrow(a_tests_db, row_num, title_str, struct('putLabels', 1))
    row_plot = set(row_plot, 'axis_labels', {'', ['row ' num2str(row_num) ]})
    plots{num_rows - row_num + 1} = row_plot;
    %border = getfield(get(row_plot, 'props'), 'border')
  else
    row_plot = plotrow(a_tests_db, row_num, title_str); %, struct('border', border)
    row_plot = set(row_plot, 'axis_labels', {'', ['row ' num2str(row_num) ]});
    plots{num_rows - row_num + 1} = row_plot;
  end
end

props(1).xLabelsPos = 'bottom';
%props.xTicksPos = 'bottom';
props.titlesPos = 'none';

tests_props = get(a_tests_db, 'props');
if isfield(tests_props, 'quiet') && tests_props.quiet == 1
  title_str = 'Rows';
else 
  title_str = get(a_tests_db, 'id');
end

a_plot = plot_stack(plots, axis_limits, orientation, title_str, props);
