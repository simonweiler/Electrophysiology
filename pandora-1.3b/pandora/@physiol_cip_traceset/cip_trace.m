function a_cip_trace = cip_trace(traceset, trace_index, props)

% cip_trace - Loads a cip_trace object from a raw data file in the traceset.
%
% Usage:
% a_cip_trace = cip_trace(traceset, trace_index, props)
%
% Description:
%
%   Parameters:
%	traceset: A physiol_cip_traceset object.
%	trace_index: Index of file in traceset.
%	props: A structure with any optional properties.
%	  showParamsList: Cell array of params to add to id field.
%	  showName: Show the name of the cell in the id field (default=1).
%	  TracesetIndex: Indicates in the id field.
%		
%   Returns:
%	a_cip_trace: A cip_trace object that holds the raw data.
%
% See also: itemResultsRow, params_tests_fileset, paramNames, testNames
%
% $Id: cip_trace.m 1335 2012-04-19 18:04:32Z cengique $
%
% Author: Cengiz Gunay <cgunay@emory.edu>, 2005/07/13
%
% Modified by:
%	Li, Su <su.li@emory.edu>, 2007/06/10 for adding HDF5 compatibility.

% Copyright (c) 2007 Cengiz Gunay <cengique@users.sf.net>.
% This work is licensed under the Academic Free License ("AFL")
% v. 3.0. To view a copy of this license, please look at the COPYING
% file distributed with this software or visit
% http://opensource.org/licenses/afl-3.0.php.

if ~ exist('props', 'var')
  props = struct;
end

% First call CIPform and get necessary parameters
[type, on, off, finish, bias, pulse] = CIPform(traceset, trace_index);
pulse_time_start=on;
pulse_time_width= off - on + 1;
ygain = 1 / traceset.vgain;
new_props=struct('type', type, 'on', on, 'off', off, 'finish', finish, 'bias', ...
		 bias, 'pulse', pulse, 'channel', traceset.vchan, 'scale_y', ...
		 ygain, 'traces', num2str(getItem(traceset, trace_index)));

if isfield(props, 'TracesetIndex')
  trace_id = [ 's' num2str(props.TracesetIndex) ',t' num2str(trace_index) ];
else
  trace_id = [ 't' num2str(trace_index) ];
end

if isfield(props, 'showParamsList')
  if ischar(props.showParamsList)
    num_params = 1;
  else
    num_params = length(props.showParamsList);
  end
  param_names = paramNames(traceset);
  param_vals = getItemParams(traceset, trace_index, ...
			     results_profile(struct, '', new_props));
  
  for param_num=1:num_params
    param_idx = strmatch(props.showParamsList{param_num}, param_names);
    trace_id = [trace_id ',' param_names{param_idx} '=' num2str(param_vals(param_idx))];
  end
end

% Put id of traceset as id
if ~ isfield(props, 'showName') || props.showName == 1
  trace_id = [get(traceset, 'neuron_id') '(' trace_id ')'];
end

traceset_props = get(traceset, 'props');
if isfield(traceset_props, 'Trials')
  % if given, take per-trial props
  trace_props = traceset_props.Trials{trace_index};
else
  % o/w use the traceset-wide props
  trace_props = traceset_props;
end
    
a_cip_trace = cip_trace(traceset.data_src, get(traceset, 'dt'), ...
			get(traceset, 'dy'), ...
			pulse_time_start, pulse_time_width, ...
			trace_id, ...
			mergeStructs(trace_props, new_props));
