function a_trace_profile = trace_profile(fileset, file_index)

% trace_profile - Loads a raw trace_profile given a file_index to this fileset.
%
% Usage:
% a_trace_profile = trace_profile(fileset, file_index)
%
% Description:
%
%   Parameters:
%	fileset: A params_tests_fileset.
%	file_index: Index of file in fileset.
%		
%   Returns:
%	a_trace_profile: A trace_profile object.
%
% See also: trace_profile, params_tests_fileset
%
% $Id: trace_profile.m 1335 2012-04-19 18:04:32Z cengique $
%
% Author: Cengiz Gunay <cgunay@emory.edu>, 2004/09/13

% Copyright (c) 2007 Cengiz Gunay <cengique@users.sf.net>.
% This work is licensed under the Academic Free License ("AFL")
% v. 3.0. To view a copy of this license, please look at the COPYING
% file distributed with this software or visit
% http://opensource.org/licenses/afl-3.0.php.

filename = getItem(fileset, file_index);
fullname = fullfile(fileset.path, filename);

% Load a trace_profile object
a_trace_profile = trace_profile(fullname, get(fileset, 'dt'), get(fileset, 'dy'), ...
				[get(fileset, 'id') '(' num2str(file_index) ')'], ...
				get(fileset, 'props'));
