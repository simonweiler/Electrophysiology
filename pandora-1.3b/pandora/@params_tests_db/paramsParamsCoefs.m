function pp_coefs = paramsParamsCoefs(a_db, p_t3ds, p_coefs)

% paramsParamsCoefs - Calculates a corrcoefs_db for each param from correlations of variant params and invariant param coefs and collects them in a cell array.
%
% Usage:
% pp_coefs = paramsParamsCoefs(a_db, p_t3ds, p_coefs)
%
% Description:
%   Skips the 'ItemIndex' test.
%
%   Parameters:
%	a_db: A tests_db object.
%	p_t3ds: Cell array of invariant parameter databases.
%	p_coefs: Cell array of tests coefficients for each parameter.
%		
%   Returns:
%	pp_coefs: A cell array of corrcoefs_dbs for each param 
%		  combination in a_db.
%
% See also: params_tests_profile
%
% $Id: paramsParamsCoefs.m 1335 2012-04-19 18:04:32Z cengique $
%
% Author: Cengiz Gunay <cgunay@emory.edu>, 2004/10/17

% Copyright (c) 2007 Cengiz Gunay <cengique@users.sf.net>.
% This work is licensed under the Academic Free License ("AFL")
% v. 3.0. To view a copy of this license, please look at the COPYING
% file distributed with this software or visit
% http://opensource.org/licenses/afl-3.0.php.

num_params = a_db.num_params;
num_tests = dbsize(a_db, 2) - num_params - 1; % Except the file indices

pp_coefs = cell([num_params, num_params]);
for param_num=1:num_params
  a_t3d = p_t3ds(param_num);
  a_coef_db = p_coefs{param_num}
  % Remove the Index columns from this
  all_test_cols = [];
  all_test_cols(1:dbsize(a_coef_db, 2)) = true(1);
  all_test_cols(tests2cols(a_coef_db, {'PageIndex'})) = false(1,1)
  joint_db = joinPages(onlyRowsTests(a_t3d, ':', 'RowIndex'), ...
                       onlyRowsTests(a_coef_db, all_test_cols, ':'));
  % Remove the Index columns from this
  all_test_cols = [];
  all_test_cols(1:dbsize(joint_db, 2)) = true(1);
  all_test_cols(tests2cols(joint_db, {'RowIndex'})) = false(1,1);
  j2_db = joinRows(onlyRowsTests(a_db, ':', 1:num_params), ...
                   onlyRowsTests(joint_db, ':', all_test_cols));
  for cparam_num=1:num_params
    all_test_cols = [];
    all_test_cols((num_params + 1):dbsize(j2_db, 2)) = true(1);
    pp_coefs{param_num, cparam_num} = ...
	corrCoefs(j2_db, param_num, find(all_test_cols), struct('skipCoefs', 0));
  end
end
