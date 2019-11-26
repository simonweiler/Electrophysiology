function a_db = unique(a_db)

% unique - Returns DB with unique rows.
%
% Usage:
% a_db = unique(a_db)
%
% Parameters:
%   a_db: tests_db from which to find uniques.
%		
% Returns:
%   a_db: The resulting tests_db.
%
% Description:
%   Keeps the original DB order.
%
% See also: unique
%
% $Id: unique.m 1334 2012-04-19 18:02:13Z cengique $
%
% Author: Cengiz Gunay <cgunay@emory.edu>, 2007/11/19

% Copyright (c) 2007 Cengiz Gunay <cengique@users.sf.net>.
% This work is licensed under the Academic Free License ("AFL")
% v. 3.0. To view a copy of this license, please look at the COPYING
% file distributed with this software or visit
% http://opensource.org/licenses/afl-3.0.php.

% indicate in name
a_db = set(a_db, 'id', [ 'uniques of ' get(a_db, 'id') ]);

% filter contents
a_db = set(a_db, 'data', uniqueValues(get(a_db, 'data')));
