function a_db = plus(left_obj, right_obj)

% plus - Adds a DB to another or to a scalar.
%
% Usage:
% a_db = plus(left_obj, right_obj)
%
% Description:
%   If DBs have mismatching columns only the common columns will be kept.
% In any case, the resulting DB columns will be sorted in the order of the
% left-hand-side DB.
%
%   Parameters:
%	left_obj, right_obj: Operands of the addition. One must be of type tests_db
%		and the other can be a scalar or tests_db.
%		
%   Returns:
%	a_db: The resulting tests_db.
%
% See also: plus
%
% $Id: plus.m 1334 2012-04-19 18:02:13Z cengique $
%
% Author: Cengiz Gunay <cgunay@emory.edu>, 2007/12/13

% Copyright (c) 2007 Cengiz Gunay <cengique@users.sf.net>.
% This work is licensed under the Academic Free License ("AFL")
% v. 3.0. To view a copy of this license, please look at the COPYING
% file distributed with this software or visit
% http://opensource.org/licenses/afl-3.0.php.

a_db = rop(left_obj, right_obj, @plus, '+');
