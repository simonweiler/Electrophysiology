function db = fillMissingParams(db, col_names, fill_value)

% fillMissingParams - Add missing columns as params with given default fill value.
%
% Usage:
% db = fillMissingParams(db, col_names, fill_value)
%
% Parameters:
%	db: A tests_db object.
%	col_names: A cell array of param names.
%	fill_value: Value to be used for missing columns.
%		
% Returns:
%	db: The tests_db object that includes the newly filled columns.
%
% Description:
%
% See also: tests_db/fillMissingColumns, params_tests_db/addParams, params_tests_db/unionCat
%
% $Id: fillMissingParams.m 1334 2012-04-19 18:02:13Z cengique $
%
% Author: Cengiz Gunay <cgunay@emory.edu>, 2008/06/02 
% - based on code from Li Su in unionCat.

% Copyright (c) 2007-2008 Cengiz Gunay <cengique@users.sf.net>.
% This work is licensed under the Academic Free License ("AFL")
% v. 3.0. To view a copy of this license, please look at the COPYING
% file distributed with this software or visit
% http://opensource.org/licenses/afl-3.0.php.

all_col_names = getParamNames(db);

% add columns in col_names but not in db.
dif_names = setdiff(col_names, all_col_names);
if ~isempty(dif_names)
    db = ...
        addColumns(db, dif_names, repmat(fill_value, dbsize(db,1), length(dif_names)));
end
