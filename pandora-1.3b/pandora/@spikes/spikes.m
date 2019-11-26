function obj = spikes(times, num_samples, dt, id)

% spikes - Spike times from a trace.
%
% Usage:
% obj = spikes(times, num_samples, dt, id)
%
%   Parameters:
%	times: The spike times [dt].
%	num_samples: Number of samples in the original trace.
%	dt: Time resolution [s].
%	id: Identification string.
%
% Description:
%		
%   Returns a structure object with the following fields:
%	times, num_samples, dt, id.
%
% General methods of spikes objects:
%   spikes		- Construct a new spikes object.
%   plot		- Graph the spikes.
%   display		- Returns and displays the identification string.
%   subsref		- Allows usage of . operator.
%   spikeRate		- Average spike rate [Hz].
%   spikeRateISI	- Average spike rate, calculated from ISIs [Hz].
%   spikeISIs		- Vector of spike ISIs [dt].
%   ISICV		- ISI coefficient of variation.
%   SFA			- Spike frequency accommodation.
%   spikeAmpSlope	- Approximate spike amplitude decay parameters.
%   withinPeriod 	- Returns a spikes object valid only within the 
%			given period.
%   periodWhole		- Period object covering the whole period.
%   getResults		- Calculates a set of tests.
%
% Additional methods:
%   See methods('spikes')
%
% See also: trace/spikes, trace, spike_shape, period
%
% $Id: spikes.m 1335 2012-04-19 18:04:32Z cengique $
%
% Author: 
%   Cengiz Gunay <cgunay@emory.edu>, 2004/07/30
%   Inspired by cip_dataset of Jeremy Edgerton.

% Copyright (c) 2007 Cengiz Gunay <cengique@users.sf.net>.
% This work is licensed under the Academic Free License ("AFL")
% v. 3.0. To view a copy of this license, please look at the COPYING
% file distributed with this software or visit
% http://opensource.org/licenses/afl-3.0.php.

if nargin == 0 % Called with no params
   obj.times = [];
   obj.num_samples = 0;
   obj.dt = 1;
   obj.id = '';
   obj = class(obj, 'spikes');
 elseif isa(times,'spikes') % copy constructor?
   obj = times;
 else
   obj.times = times;
   obj.num_samples = num_samples;
   obj.dt = dt;
   obj.id = id;
   obj = class(obj, 'spikes');
end

