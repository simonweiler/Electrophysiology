function p_t3ds = invarParams(a_db)

% invarParams - Calculates invariant param dbs for all parameters and returns in an array.
%
% Usage:
% p_t3ds = invarParams(a_db)
%
% Description:
%   Skips the 'ItemIndex' test.
%
%   Parameters:
%	a_db: A tests_db object.
%		
%   Returns:
%	p_t3ds: An array of tests_3D_dbs for each param in a_db.
%
% See also: params_tests_profile
%
% $Id: invarParams.m 1335 2012-04-19 18:04:32Z cengique $
%
% Author: Cengiz Gunay <cgunay@emory.edu>, 2004/10/17

% Copyright (c) 2007 Cengiz Gunay <cengique@users.sf.net>.
% This work is licensed under the Academic Free License ("AFL")
% v. 3.0. To view a copy of this license, please look at the COPYING
% file distributed with this software or visit
% http://opensource.org/licenses/afl-3.0.php.

num_params = a_db.num_params;
num_tests = dbsize(a_db, 2) - num_params; 

% no need, already done in invarValues
%a_db = sortrows(a_db, 1:num_params);

p_t3ds(1:num_params) = tests_3D_db;
for param_num=1:num_params
  p_t3ds(param_num) = invarParam(a_db, param_num);
end
