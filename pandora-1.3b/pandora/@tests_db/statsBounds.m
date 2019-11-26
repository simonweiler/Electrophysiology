function a_stats_db = statsBounds(a_db, tests, props)

% statsBounds - Generates a stats_db object with three rows corresponding to the mean, min, max and number of observations of the tests' distributions. 
%
% Usage:
% a_stats_db = statsBounds(a_db, tests, props)
%
% Description:
%   A page is generated for each page of data in db.
%
%   Parameters:
%	a_db: A tests_db object.
%	tests: A selection of tests (see onlyRowsTests).
%	props: A structure with any optional properties for stats_db.
%		
%   Returns:
%	a_stats_db: A stats_db object.
%
% See also: tests_db
%
% $Id: statsBounds.m 1335 2012-04-19 18:04:32Z cengique $
%
% Author: Cengiz Gunay <cgunay@emory.edu>, 2004/10/07

% Copyright (c) 2007 Cengiz Gunay <cengique@users.sf.net>.
% This work is licensed under the Academic Free License ("AFL")
% v. 3.0. To view a copy of this license, please look at the COPYING
% file distributed with this software or visit
% http://opensource.org/licenses/afl-3.0.php.

if ~ exist('props', 'var')
  props = struct([]);
end

if ~ exist('tests', 'var')
  tests = ':';
end

cols = tests2cols(a_db, tests);

num_pages = dbsize(a_db, 3);
pages=1:num_pages;
data = repmat(NaN, [4, length(cols), num_pages]);
for page_num=pages
  a_page_db = onlyRowsTests(a_db, ':', tests, page_num);
  if dbsize(a_page_db, 1) > 0
    [means, n] = mean(a_page_db, 1);
    data(:, :, page_num) = [get(means, 'data'); ...
                        min(get(a_page_db, 'data'), [], 1); ...
                        max(get(a_page_db, 'data'), [], 1); n];
  end
end

row_names = {'mean', 'min', 'max', 'n'};

% Original column names
col_name_cell = fieldnames(a_db.col_idx);
col_names = col_name_cell(cols);

a_stats_db = stats_db(data, col_names, row_names, {}, ...
		      [ 'Mean value and min-max bounds from ' get(a_db, 'id') ], props);
