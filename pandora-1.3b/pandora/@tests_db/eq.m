function rows = eq(db, row)

% eq - Equality (==) operator. Returns logical indices of db rows 
%	that match with given row.
%
% Usage:
% rows = eq(db, row)
%
% Description:
%
%   Parameters:
%	db: A tests_db object.
%	row: Row array to be compared with db rows.
%		
%   Returns:
%	rows: A logical or index vector to be used in indexing db objects. 
%
% See also: eq, tests_db
%
% $Id: eq.m 1335 2012-04-19 18:04:32Z cengique $
%
% Author: Cengiz Gunay <cgunay@emory.edu>, 2004/09/17

% Copyright (c) 2007 Cengiz Gunay <cengique@users.sf.net>.
% This work is licensed under the Academic Free License ("AFL")
% v. 3.0. To view a copy of this license, please look at the COPYING
% file distributed with this software or visit
% http://opensource.org/licenses/afl-3.0.php.

% - find 0 rows in db
rows = ~ compareRows(db, row);
