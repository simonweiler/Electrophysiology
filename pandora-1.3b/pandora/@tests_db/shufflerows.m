function a_db = shufflerows(db, tests, grouped)

% shufflerows - Returns a db with rows of given test columns are shuffled. 
%
% Usage:
% s = shufflerows(db, tests, grouped)
%
% Description:
%   Can be used for shuffle prediction. Basically, shuffle rows of tests and rerun
% high order analyses. 
%
%   Parameters:
%	db: A tests_db object.
%	tests: Tests to shuffle.
%	grouped: If 1 then shuffle tests all together, 
%		if 0 shuffle each test separately (default=0).
%		
%   Returns:
%	a_db: The shuffled db.
%
% See also: tests_db
%
% $Id: shufflerows.m 1335 2012-04-19 18:04:32Z cengique $
%
% Author: Cengiz Gunay <cgunay@emory.edu>, 2004/11/10

% Copyright (c) 2007 Cengiz Gunay <cengique@users.sf.net>.
% This work is licensed under the Academic Free License ("AFL")
% v. 3.0. To view a copy of this license, please look at the COPYING
% file distributed with this software or visit
% http://opensource.org/licenses/afl-3.0.php.

if ~ exist('grouped', 'var')
  grouped = 0;
end

cols = tests2cols(db, tests);
num_rows = dbsize(db, 1);
data = get(db, 'data');

if grouped == 1
    data(:, cols, 1) = data(randperm(num_rows)', cols, 1);
else
  for col_num = cols
    data(:, col_num, 1) = data(randperm(num_rows)', col_num, 1);
  end
end

a_db = set(db, 'data', data);