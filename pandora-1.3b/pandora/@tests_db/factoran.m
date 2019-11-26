function a_factors_db = factoran(db, num_factors, props)

% factoran - Generates a database of factor loadings obtained from the 
%		factor analysis of db with factoran. Each row corresponds
%		to a rotated factor and columns represent observed variables.
%
% Usage:
% a_factors_db = factoran(db, num_factors, props)
%
% Description:
%  Uses the promax method to rotate common factors.
%
%   Parameters:
%	db: A tests_db object.
%	num_factors: Number of common factors to look for.
%	props: A structure with any optional properties.
%		
%   Returns:
%	a_factors_db: A corrcoefs_db of the coefficients and page indices.
%
% See also: tests_db, corrcoefs_db
%
% $Id: factoran.m 1335 2012-04-19 18:04:32Z cengique $
%
% Author: Cengiz Gunay <cgunay@emory.edu>, 2004/11/08

% Copyright (c) 2007 Cengiz Gunay <cengique@users.sf.net>.
% This work is licensed under the Academic Free License ("AFL")
% v. 3.0. To view a copy of this license, please look at the COPYING
% file distributed with this software or visit
% http://opensource.org/licenses/afl-3.0.php.

if ~ exist('props', 'var')
  props = struct([]);
end

[LoadingsPM, specVarPM, T, stats] = ...
    factoran(get(db, 'data'), num_factors, 'rotate', 'promax' );


a_factors_db = tests_db(LoadingsPM', fieldnames(get(db, 'col_idx')), {}, ...
			[ 'Common factors of ' get(db, 'id') ], props);

