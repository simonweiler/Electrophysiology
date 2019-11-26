function [avg_tr sd_tr] = avgTraces(traces, props)

% avgTraces - Average multiple traces.
%
% Usage: 
% [avg_tr sd_tr] = avgTraces(traces, props)
%
% Parameters:
%   traces: A vector of trace objects.
%   props: A structure with any optional properties.
%     calcSE: If given, calculate standard error instead of deviation.
%     id: String to replace the id property of averaged trace. The term
%         "average" or "SD" will be prepended to it.  By default
%         it will be lengthy and show the arithmetic done.
%
% Returns:
%	avg_tr: A trace object that holds the average.
%	sd_tr: A trace object that holds the standard deviation or error.
%
% Description:
%
% See also: trace
%
% $Id: avgTraces.m 1334 2012-04-19 18:02:13Z cengique $
%
% Author: Cengiz Gunay <cgunay@emory.edu>, 2010/11/09

% Copyright (c) 2010 Cengiz Gunay <cengique@users.sf.net>.
% This work is licensed under the Academic Free License ("AFL")
% v. 3.0. To view a copy of this license, please look at the COPYING
% file distributed with this software or visit
% http://opensource.org/licenses/afl-3.0.php.

% manually tested, make a unit test in the future

if ~ exist('props', 'var')
  props = struct;
end

% TODO: implement separately as the sum() function
num_traces = length(traces);
if num_traces > 0
  avg_tr = traces(1);  
  for trace_num = 2:num_traces
      avg_tr = avg_tr + traces(trace_num);
  end
  avg_tr = avg_tr ./ num_traces;
  if isfield(props, 'id')
    avg_tr = set(avg_tr, 'id', ['Avg-' props.id ]);
  end
end

% SD/E
if num_traces > 0
  sd_tr = (traces(1) - avg_tr).^2;
  for trace_num = 2:num_traces
      sd_tr = sd_tr + (traces(trace_num) - avg_tr).^2;
  end
  sd_tr = sd_tr ./ num_traces;
  if isfield(props, 'calcSE')
      sd_tr = sd_tr ./ num_traces;
  end
  sd_tr = sqrt(sd_tr);
  if isfield(props, 'id')
    sd_tr = set(sd_tr, 'id', ['SD-' props.id ]);
  end
end