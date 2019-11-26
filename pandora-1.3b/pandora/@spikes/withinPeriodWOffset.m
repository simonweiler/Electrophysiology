function obj = withinPeriodWOffset(s, a_period)

% withinPeriodWOffset - Returns a spikes object valid only within the given period, keeps the offset.
%
% Usage:
% obj = withinPeriodWOffset(s, a_period)
%
% Description:
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
% $Id: withinPeriodWOffset.m 1335 2012-04-19 18:04:32Z cengique $
%
% Author: Cengiz Gunay <cgunay@emory.edu>, 2005/05/09

% Copyright (c) 2007 Cengiz Gunay <cengique@users.sf.net>.
% This work is licensed under the Academic Free License ("AFL")
% v. 3.0. To view a copy of this license, please look at the COPYING
% file distributed with this software or visit
% http://opensource.org/licenses/afl-3.0.php.

% check if input is an array
num_objs = length(s);
if num_objs > 1
  obj = repmat(spikes, 1, num_objs);
  for obj_num = 1:num_objs
    obj(obj_num) = ...
        withinPeriodWOffset(s(obj_num), a_period);
  end
  return
end

% for single spikes object
s.times = s.times(s.times > a_period.start_time & s.times <= a_period.end_time);
s.num_samples = a_period.end_time - a_period.start_time;

obj = s;

