function param_names = paramNames(fileset)

% paramNames - Returns the only parameter, 'pAcip,' for this fileset.
%
% Usage:
% param_names = paramNames(fileset)
%
% Description:
% Looks at the filename of the first file to find the parameter names.
%
%   Parameters:
%	fileset: A params_tests_fileset.
%		
%   Returns:
%	params_names: Cell array with ordered parameter names.
%
% See also: params_tests_fileset, paramNames, testNames
%
% $Id: paramNames.m 1335 2012-04-19 18:04:32Z cengique $
%
% Author: Cengiz Gunay <cgunay@emory.edu>, 2004/12/06

% Copyright (c) 2007 Cengiz Gunay <cengique@users.sf.net>.
% This work is licensed under the Academic Free License ("AFL")
% v. 3.0. To view a copy of this license, please look at the COPYING
% file distributed with this software or visit
% http://opensource.org/licenses/afl-3.0.php.

% CIP magnitude in pA
param_names = { 'NeuronId', 'pAcip' };