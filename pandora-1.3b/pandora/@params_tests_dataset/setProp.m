function obj = setProp(obj, varargin)
% setProp - Generic method for setting optional object properties.
%
% Usage:
% obj = setProp(obj, prop1, val1, prop2, val2, ...)
%
% Description:
%   Modifies or adds property values. As many property name-value 
% pairs can be specified.
%
%   Parameters:
%	obj: Any object that has a props field.
%	attr: Property name
%	val: Property value.
%
%   Returns:
%	obj: The new object with the updated properties.
%
% See also: 
%
% $Id: setProp.m 1334 2012-04-19 18:02:13Z cengique $
%
% Author: Cengiz Gunay <cgunay@emory.edu>, 2004/11/22

% Copyright (c) 2007 Cengiz Gunay <cengique@users.sf.net>.
% This work is licensed under the Academic Free License ("AFL")
% v. 3.0. To view a copy of this license, please look at the COPYING
% file distributed with this software or visit
% http://opensource.org/licenses/afl-3.0.php.

for index=1:2:length(varargin)
  attr = varargin{index};
  val = varargin{index + 1};
  try
    a = obj.props;
    a(1).(attr) = val;
    obj = set(obj, 'props', a);
  catch
    rethrow(lasterror);
  end
end
