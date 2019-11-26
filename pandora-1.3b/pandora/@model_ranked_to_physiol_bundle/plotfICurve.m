function a_doc = docfICurve(r_bundle, rank_num, props)

% docfICurve - OBSOLETE - Generates a a_p of given rank_num from the ranked_bundle.
%
% Usage:
% a_doc = docfICurve(r_bundle, crit_bundle, crit_db, props)
%
% Description:
%
%   Parameters:
%	r_bundle: A ranked_bundle object.
%	rank_num: Rank index for which to generate the a_doc.
%	props: A structure with any optional properties.
%
%   Returns:
%	a_doc: A doc_plot that contains a f-I curve plot and associated captions.
%
%   Example:
% >> a_d = docfICurve(r, 1);
% >> plot(a_d, 'The f-I curve of best matching model');
%
% See also: doc_generate, doc_plot
%
% $Id: plotfICurve.m 1335 2012-04-19 18:04:32Z cengique $
%
% Author: Cengiz Gunay <cgunay@emory.edu>, 2006/01/16

% Copyright (c) 2007 Cengiz Gunay <cengique@users.sf.net>.
% This work is licensed under the Academic Free License ("AFL")
% v. 3.0. To view a copy of this license, please look at the COPYING
% file distributed with this software or visit
% http://opensource.org/licenses/afl-3.0.php.

a_ranked_db = r_bundle.ranked_db;
joined_db = joinOriginal(a_ranked_db, 1); % Get only the best

% adjust labels for LaTeX
a_db_id = lower(properTeXLabel(get(a_ranked_db.orig_db, 'id')));
crit_db_id = lower(properTeXLabel(get(a_ranked_db.crit_db, 'id')));
traceset_index = a_ranked_db.crit_db(1, 'TracesetIndex', 1).data;
crit_trace_id = lower(properTeXLabel([ get(getItem(r_bundle.crit_bundle.dataset, traceset_index), 'id') ...
				      '(s' num2str(traceset_index) ')']));

short_caption = [ 'f-I curves of best matching model to ' crit_db_id '.' ];
caption = [ short_caption ];
curve_pAvals = [0 40 100 200];
curve_tests = {'IniSpontSpikeRateISI_0pA', 'PulseIni100msSpikeRateISI_D40pA', ...
	       'PulseIni100msSpikeRateISI_D100pA', 'PulseIni100msSpikeRateISI_D200pA'};
curve_labels = {'current pulse [pA]', 'firing rate [Hz]'};
best_trial_num = joined_db(1, 'trial').data;
a_doc = ...
    doc_plot(plot_superpose({plotYTests(statsMeanStd(r_bundle.crit_bundle.joined_db), ...
					curve_pAvals, curve_tests, curve_labels, ...
					'', 'phys. avg.', [], struct('quiet', 1)), ...
			     plotYTests(a_ranked_db.crit_db(1, :), ...
					curve_pAvals, curve_tests, ...
					curve_labels, '', ...
					[ crit_trace_id ' (avg)'], ...
					[], struct('quiet', 1)), ...
			     plotYTests(joined_db(1, :), curve_pAvals, curve_tests, ...
					curve_labels, '', ...
					[ 'model (t' num2str(best_trial_num) ')'], ...
					[], struct('quiet', 1))}, {}, ...
			    'f-I curve of best matching model'), ...
	     caption, [crit_db_id  ' - fI curve of best matching model from ' a_db_id ], ...
	     struct('floatType', 'figure', 'center', 1, ...
		    'width', '.7\textwidth', 'shortCaption', short_caption), ...
	     'frequency-current curve', struct);
