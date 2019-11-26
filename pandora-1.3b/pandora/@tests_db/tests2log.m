function a_log = tests2log(db, dim_num, tests)

% tests2log - Return logical array of indices from a test names/numbers specification.
%
% Usage:
% a_log = tests2log(db, dim_name, tests)
%
% Parameters:
%	db: A tests_db object.
%	dim_num: Number between 1-3 to choose dimension: row, column, or page.
%	tests: Either a single or array of column numbers, or a single
%		test name or a cell array of test names. If ':', all
%		tests. For name strings, regular expressions are
%		supported if quoted with slashes (e.g., '/a.*/')
%		
% Returns:
%	a_log: Array of column indices.
%
% Description:
%   See tests2idx.
%
% Example:
% >> cols = tests2log(a_db, 'col', {'col1', '/col2+/'});
% >> stripped_db = a_db(:, ~cols)
% will remove columns col1 and col2, col22, col22, etc. from stripped_db.
%
% See also: tests_db, tests2cols, regexp
%
% $Id: tests2log.m 1334 2012-04-19 18:02:13Z cengique $
%
% Author: Cengiz Gunay <cgunay@emory.edu>, 2008/05/27

% Copyright (c) 2007 Cengiz Gunay <cengique@users.sf.net>.
% This work is licensed under the Academic Free License ("AFL")
% v. 3.0. To view a copy of this license, please look at the COPYING
% file distributed with this software or visit
% http://opensource.org/licenses/afl-3.0.php.

% get logical indices from indices
a_log = false(dbsize(db, dim_num), 1);
a_log(tests2idx(db, dim_num, tests)) = true;