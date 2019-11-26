function obj = withinPeriod(s, a_period)

% withinPeriod - Returns a spikes object valid only within the given period, subtracts the offset.
%
% Usage:
% obj = withinPeriod(s, a_period)
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
% $Id: withinPeriod.m 1335 2012-04-19 18:04:32Z cengique $
%
% Author: Cengiz Gunay <cgunay@emory.edu>, 2004/07/31
% Modified:

% Copyright (c) 2007 Cengiz Gunay <cengique@users.sf.net>.
% This work is licensed under the Academic Free License ("AFL")
% v. 3.0. To view a copy of this license, please look at the COPYING
% file distributed with this software or visit
% http://opensource.org/licenses/afl-3.0.php.

% TODO:
% - Relate this method by overloading an abstract class/interface periodable(?) 

s = withinPeriodWOffset(s, a_period);

% reset the offset
% check if input is an array
num_objs = length(s);
if num_objs > 1
  for obj_num = 1:num_objs
    if ~ isempty(s(obj_num).times)
      s(obj_num).times = s(obj_num).times - a_period.start_time + 1; 
    end
  end
else
  if ~ isempty(s.times)
    s.times = s.times - a_period.start_time + 1; 
  end
end

obj = s;

