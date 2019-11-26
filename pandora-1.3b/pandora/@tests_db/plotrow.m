function a_plot = plotrow(a_tests_db, row, title_str, props)

% plotrow - Creates a plot_abstract describing the desired db row.
%
% Usage:
% a_plot = plotrow(a_tests_db, row, title_str, props)
%
% Description:
%
%   Parameters:
%	a_tests_db: A tests_db object.
%	row: Row number to visualize.
%	title_str: Optional title string.
%	props: A structure with any optional properties.
%	  putLabels: Put special column name labels.
%		
%   Returns:
%	a_plot: A plot_abstract object that can be plotted.
%
% See also: plot_abstract, plotFigure
%
% $Id: plotrow.m 1335 2012-04-19 18:04:32Z cengique $
%
% Author: Cengiz Gunay <cgunay@emory.edu>, 2004/11/08

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

data = a_tests_db.data;
x_vals = 1:dbsize(a_tests_db, 2);
%props.XTickLabel = fieldnames(get(a_tests_db, 'col_idx'));
props.XTick = x_vals;
props.grid = 1;

if isfield(props, 'putLabels') && props.putLabels == 1
  props.XTickLabel = {''};
end

% TODO: need special plot_abstract for making colored bars:
%rows = [x_vals, data(row, :, 1)];
%flatrow = num2cell(reshape(rows, 1, 2*length(x_vals)));

if isfield(props, 'quiet')
  the_title = title_str;
else
  the_title = [ get(a_tests_db, 'id') title_str ];
end

% test a special case where bar fails
if length(row) > 1 && size(data, 2) == 1
  x_vals = [ x_vals; (length(x_vals) + 1) ];
  data = [ data, repmat(NaN, length(row), 1) ];
end

% get row names if available
row_idx = get(a_tests_db, 'row_idx');
if ~ isempty(fieldnames(row_idx))
  all_row_names = fieldnames(row_idx);
  row_names = [ all_row_names(row) ];
else
  row_names = {};
end

a_plot = ...
    plot_abstract({ x_vals, data(row, :, 1)' }, {'', ''}, ...
		  properTeXLabel(the_title), properTeXLabel(row_names), 'bar', ...
		  mergeStructs(props, struct('tightLimits', 1)));

add_props = struct;
if isfield(props, 'putLabels') && props.putLabels == 1
  label_plots = cell(1, dbsize(a_tests_db, 2));

  % align label with axis limits if given
  if isfield(props, 'axisLimits') && ...
      ~isnan(props.axisLimits(3)) && ~isinf(props.axisLimits(3))
    min_val = props.axisLimits(3);
  else
    % Bars start from zero if no
    % negative bars exist
    min_val = min([min(data(row, :, 1)), 0]);  
  end

  % to make space for labels
  add_props = ...
        struct('border', [0.05 0.2 0 0]);
  
  col_names = fieldnames(get(a_tests_db, 'col_idx'));
  for col=1:dbsize(a_tests_db, 2)
    label_plots{col} = ...
	plot_abstract({ col, min_val, properTeXLabel(col_names{col}), ...
		       struct('HorizontalAlignment', 'right', ...
                              'VerticalAlignment', 'top', ...
			      'Rotation', 20)}, {'', ''}, ...
		      properTeXLabel(get(a_tests_db, 'id')), {}, ...
                      'subTextLabel', props);
% $$$     plot_abstract({ col/(dbsize(a_tests_db, 2) + 1), -0.03, properTeXLabel(col_names{col}), ...
% $$$ 		       struct('Units', 'normalized', ...
% $$$ 			      'HorizontalAlignment', 'right', ...
% $$$ 			      'Rotation', 20)}, {'', ''}, ...
% $$$ 		      properTeXLabel(get(a_tests_db, 'id')), {}, 'subTextLabel', props);
  end
  a_plot = plot_superpose( [{a_plot}, label_plots], {}, '', ...
                           mergeStructs(props, add_props));
end