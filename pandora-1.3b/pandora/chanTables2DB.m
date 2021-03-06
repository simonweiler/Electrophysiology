function a_chans_db = chanTables2DB(tables, id, props)

% chanTables2DB - Creates a DB with channel tables exported from Genesis.
%
% Usage: 
% a_chans_db = chanTables2DB(tables, id, props)
%
% Description:
%
%   Parameters:
%	tables: Structures returned from the dump files generated by dump_chans.g.
%	id: String that identify the source of the tables structure.
%	props: A structure with any optional properties.
%	  (rest passed to tests_db.)
%
%   Returns:
%	a_chans_db: A tests_db object containing channel tables.
%
% See also: trace, trace/plot, plot_abstract, GP/common/dump_chans.g (Genesis)
%
% $Id: chanTables2DB.m 1335 2012-04-19 18:04:32Z cengique $
%
% Author: Cengiz Gunay <cgunay@emory.edu>, 2007/03/07

% Copyright (c) 2007 Cengiz Gunay <cengique@users.sf.net>.
% This work is licensed under the Academic Free License ("AFL")
% v. 3.0. To view a copy of this license, please look at the COPYING
% file distributed with this software or visit
% http://opensource.org/licenses/afl-3.0.php.

if ~ exist('props', 'var')
  props = struct;
end

chan_names = fieldnames(tables)';
a_db = tests_db;
channel_info = struct;

chan_num = 1;
% go thru all channels in tables
for chan_name = chan_names
  chan_name = chan_name{1};

  chan = tables.(chan_name);
  chan_fields = fieldnames(chan)';
  gate_names = chan_fields(~cellfun(@isempty, regexp(chan_fields, '.*_minf|.*_tau', 'match')));
  
  %   create DB object & concat
  a_db = addColumns(a_db, makeChanDB);

  %   separate plot for each gate
  chan_num = chan_num + 1;
end

props.chan_names = chan_names;

% set the props at the end
a_chans_db = chans_db(a_db, {}, channel_info, id, props);

% inner function: return all gates of one channel in a tests_db object
function a_db = makeChanDB
  results = repmat(NaN, size(chan.(gate_names{1}), 1), length(gate_names) + 1);

  % the x-axis is always the same
  a_result = chan.(gate_names{1});
  results(:, 1) = a_result(:, 1);

  gate_num = 2;
  for gate_name = gate_names
    gate_name = gate_name{1};

    a_result = chan.(gate_name);

    results(:, gate_num) = a_result(:, 2);

    new_gate_names{gate_num - 1} = [ chan_name '_' gate_name ];
    
    gate_num = gate_num + 1;
  end

  % create the chan DB
  a_db = tests_db(results, { [ chan_name '_x' ], new_gate_names{:} }, {}, ...
                  [ id ', ' chan_name ]);

  % other fields (such as Gbar and powers) go to channel_info (assuming they are scalars)
  other_fields = setdiff(chan_fields, gate_names);
  for field_name = other_fields
    channel_info.([chan_name '_' field_name{1}]) = chan.(field_name{1});
  end


end

end

