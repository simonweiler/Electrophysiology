function obj = period(start_time, end_time)

% period - Start and end times of a period in terms of the dt of the trace to which belongs.
%
% Usage:
% obj = period(start_time, end_time)
%
% Parameters:
%   start_time, end_time: Inclusive period [dt]. If end_time is missing
%   	and start_time has two values, the second one is used as end_time.
%		
% Returns a structure object with the following fields:
%   start_time, end_time.
%
% Description:
%
% General operations on period objects:
%   period		- Construct a new period object.
%   display		- Returns and displays the identification string.
%   get			- Gets attributes of this object and parents.
%   subsref		- Allows usage of . operator.
%
% Additional methods:
%	See methods('period')
%
% See also: trace, spikes, spike_shape
%
% $Id: period.m 1335 2012-04-19 18:04:32Z cengique $
%
% Author: Cengiz Gunay <cgunay@emory.edu>, 2004/07/30

% Copyright (c) 2007 Cengiz Gunay <cengique@users.sf.net>.
% This work is licensed under the Academic Free License ("AFL")
% v. 3.0. To view a copy of this license, please look at the COPYING
% file distributed with this software or visit
% http://opensource.org/licenses/afl-3.0.php.

% TODO:
% - Maybe period should have its own dt?

if nargin == 0 % Called with no params
   obj.start_time = 0;
   obj.end_time = 0;
   obj = class(obj,'period');
 elseif isa(start_time,'period') % copy constructor?
   obj = start_time;
else
  if length(start_time) == 2 && ~ exist('end_time', 'var')
  end_time = start_time(2);
  start_time = start_time(1);
  end
  obj.start_time = start_time;
  obj.end_time = end_time;
  obj = class(obj,'period');
end