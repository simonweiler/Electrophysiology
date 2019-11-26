function a_p = plotScatter3D(a_db, test1, test2, test3, title_str, short_title, props)

% plotScatter3D - Create a 3D scatter plot of the given three tests.
%
% Usage:
% a_p = plotScatter3D(a_db, test1, test2, test3, title_str, short_title, props)
%
% Description:
%
%   Parameters:
%	a_db: A params_tests_db object.
%	test1, test2, test3: X, Y, & Z variables.
%	title_str: (Optional) A string to be concatanated to the title.
%	short_title: (Optional) Few words that may appear in legends of multiplot.
%	props: A structure with any optional properties.
%	  LineStyle: Plot line style to use. (default: 'x')
%	  Regress: Calculate and plot a linear regression.
%	  quiet: If 1, don't include database name on title.
%		
%   Returns:
%	a_p: A plot_abstract.
%
% See also: 
%
% $Id: plotScatter3D.m 1334 2012-04-19 18:02:13Z cengique $
%
% Author: Cengiz Gunay <cgunay@emory.edu>, 2007/11/30

% Copyright (c) 2007 Cengiz Gunay <cengique@users.sf.net>.
% This work is licensed under the Academic Free License ("AFL")
% v. 3.0. To view a copy of this license, please look at the COPYING
% file distributed with this software or visit
% http://opensource.org/licenses/afl-3.0.php.

if ~ exist('title_str', 'var')
  title_str = '';
end

if ~ exist('props', 'var')
  props = struct;
end

col1 = tests2cols(a_db, test1);
col2 = tests2cols(a_db, test2);
col3 = tests2cols(a_db, test3);

col1_db = onlyRowsTests(a_db, ':', col1);
col2_db = onlyRowsTests(a_db, ':', col2);
col3_db = onlyRowsTests(a_db, ':', col3);

% skip NaN value rows
non_nans_rows = ~(isnan(col1_db) | isnan(col2_db));
col1_db = onlyRowsTests(col1_db, non_nans_rows, ':');
col2_db = onlyRowsTests(col2_db, non_nans_rows, ':');
col3_db = onlyRowsTests(col3_db, non_nans_rows, ':');

test_names = fieldnames(get(a_db, 'col_idx'));

if ~ exist('short_title', 'var') || isempty(short_title)
  short_title = [strrep(test_names{col1}, '_', ' ') ', ' ...
		 strrep(test_names{col2}, '_', ' ')  ', ' ...
                 strrep(test_names{col2}, '_', ' ')];
end

if ~ isfield(props, 'quiet')
  all_title = [ strrep(get(a_db, 'id'), '_', '\_') title_str ];
else
  all_title = title_str;
end


if isfield(props, 'LineStyle')
  line_style = {props.LineStyle};
else
  line_style = {};
  props.LineStyleOrder = {'x', '+', 'd', 'o', '*', 's'};
end

if isfield(props, 'Regress')
  [b,bint,r,rint,stats] = ...
      regress(get(col3_db, 'data'), [ones(dbsize(col1_db, 1), 1), ...
                      get(col1_db, 'data') get(col2_db, 'data')]);
  if ~isempty(all_title)
    all_title = [ all_title, '; '];
  end
  all_title = [ all_title, 'regress p=' sprintf('%.4f', stats(3)) ];
end

col_labels = strrep({test_names{[col1 col2 col3]}}, '_', ' ');
a_p = plot_abstract({get(col1_db, 'data'), get(col2_db, 'data'), get(col3_db, 'data'), line_style{:}}, ...
		    { col_labels{:} }, ...
		    all_title, { short_title }, 'plot3', ...
		    props); 

if isfield(props, 'Regress')
  x_lims = [min(get(col1_db, 'data')) max(get(col1_db, 'data'))];
  y_lims = [min(get(col2_db, 'data')) max(get(col2_db, 'data'))];
  % TODO: to be completed
  a_p = plot_superpose([a_p, plot_abstract({x_lims, x_lims * b(2) + b(1), 'm-'}, ...
					   { }, '', { '' }, 'patch', props)], {}, '');
end