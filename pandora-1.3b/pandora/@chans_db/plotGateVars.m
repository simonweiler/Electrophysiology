function a_subplot = plotGateVars(a_chans_db, chan_name, gate_subnames, title_str, props)
  
% plotGateVars - Plot given channel gate variables of the same channel superposed.
%
% Usage: 
% a_plot = plotGateVars(a_chans_db, chan_name, gate_subnames, title_str, props)
%
% Description:
%
%   Parameters:
%	a_chans_db: A chans_db describing channel variables.
%	chan_name: Name of channel that make up the stem of variable
%		names.
% 	gate_subnames: Gate names of the channel.
%	title_str: (Optional) A string to be concatanated to the title.
%	props: A structure with any optional properties.
%	  usePowers: Use the gate powers, Luke.
%	  (rest passed to plot_abstract.)
%
%   Returns:
%	a_plot: A plot_abstract object that can be visualized.
%
% See also: trace, trace/plot, plot_abstract
%
% $Id: plotGateVars.m 1335 2012-04-19 18:04:32Z cengique $
%
% Author: Cengiz Gunay <cgunay@emory.edu>, 2007/07/01

% Copyright (c) 2007 Cengiz Gunay <cengique@users.sf.net>.
% This work is licensed under the Academic Free License ("AFL")
% v. 3.0. To view a copy of this license, please look at the COPYING
% file distributed with this software or visit
% http://opensource.org/licenses/afl-3.0.php.

  if ~ exist('props', 'var')
    props = struct;
  end
  
  subplots = {};
  db_id = get(a_chans_db, 'id');

  if ~exist('title_str', 'var') || isempty(title_str)
    title_str = [ strrep(chan_name, '_', ' ') ];
  end
  
  for gate_subname = gate_subnames
    % wrapped in two cells by regexp??
    gate_subname = gate_subname{1};

    if ~ isempty(gate_subname)

      if isfield(props, 'usePowers')
        power_name = [ gate_subname '_powered' ];

        % take power of channel gate
        a_chans_db = ...
            addColumns(a_chans_db, { power_name }, ...
                       get(onlyRowsTests(a_chans_db, ':', gate_subname), ...
                           'data') .^ ...
                       a_chans_db.channel_info.(strrep(gate_subname, ...
                                                       '_minf', 'power')) );
        
        gate_subname = power_name;
      end
      
      gate_subname_label = strrep(gate_subname, '_', ' ');
      subplots = { subplots{:}, ...
		  plotScatter(a_chans_db, [ chan_name '_x' ], gate_subname, title_str, ...
			      [ db_id ', ' gate_subname_label], ...
			      mergeStructs(props, struct('LineStyle', '-', 'quiet', 1)))};
    end

  end

  a_subplot = plot_superpose(subplots, {}, '');
end
