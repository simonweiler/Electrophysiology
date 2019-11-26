function obj = intoPeriod(s, a_period)

% intoPeriod - Shifts the spikes times to be within the given period.
%
% Usage:
% obj = intoPeriod(s, a_period)
%
% Description:
%   Assuming this spikes object's length fits into the given period, it shifts
% all times to start from the beginning of the given period. This may be used
% to reconstruct the original spikes object from subperiods that were cut out
% previously, using the withinPeriod method.
%
%   Parameters:
%	s: A spikes object.
%	a_period: The desired period 
%
%   Returns:
%	obj: A spikes object
%
% See also: spikes, period
%
% $Id: intoPeriod.m 1335 2012-04-19 18:04:32Z cengique $
%
% Author: Cengiz Gunay <cgunay@emory.edu>, 2004/07/31
% Modified: (see CVS log)

% Copyright (c) 2007 Cengiz Gunay <cengique@users.sf.net>.
% This work is licensed under the Academic Free License ("AFL")
% v. 3.0. To view a copy of this license, please look at the COPYING
% file distributed with this software or visit
% http://opensource.org/licenses/afl-3.0.php.

% TODO:
% - Relate this method by overloading an abstract class/interface periodable(?) 

% shift the offset
s.times = s.times + a_period.start_time - 1; 

if max(s.times) > a_period.end_time
  error('Spikes object does not fit into desired period.');
end

obj = s;

