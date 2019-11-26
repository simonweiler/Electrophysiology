function [a_db, varargout] = mean(a_db, dim)

% mean - Returns the mean of the data matrix of a_db. Ignores NaN values.
%
% Usage:
% [a_db, n] = mean(a_db, dim)
%
% Description:
%   Does a recursive operation over dimensions in order to remove NaN values.
% This takes more time than a straightforward mean operation. 
%
%   Parameters:
%	a_db: A tests_db object.
%	dim: Work down dimension.
%		
%   Returns:
%	a_db: The DB with one row of mean values.
%	n: (Optional) Numbers of non-NaN rows included in calculating each column.
%
% See also: mean, tests_db
%
% $Id: mean.m 1335 2012-04-19 18:04:32Z cengique $
%
% Author: Cengiz Gunay <cgunay@emory.edu>, 2004/10/06

% Copyright (c) 2007 Cengiz Gunay <cengique@users.sf.net>.
% This work is licensed under the Academic Free License ("AFL")
% v. 3.0. To view a copy of this license, please look at the COPYING
% file distributed with this software or visit
% http://opensource.org/licenses/afl-3.0.php.

if ~ exist('dim', 'var')
  dim = 1; % Go down rows by default
end

% Always do row-wise
order = 1:length(dbsize(a_db));
if dim ~= 1
  order(dim) = 1;
  order(1) = dim;
  data = permute(a_db.data, order);
else
  data = a_db.data;
end

% Allocate results array
db_size = size(data);
s = repmat(NaN, [1 db_size(2:end)]);

% Do a loop over EACH other dimension (!)
[s, n] = recmean(data, length(db_size));

if dim ~= 1
  s = ipermute(s, order);
end

a_db = set(a_db, 'id', [ 'Mean of ' get(a_db, 'id') ]);
a_db = set(a_db, 'data', s);

nout = max(nargout,1) - 1;

if nout > 0
  varargout{1} = n;
end

% Recursive std needed for stripping NaNs in each dimension
function [s, n] = recmean(data, dim)
  if dim == 1
    sdata = data(~isnan(data(:)) & ~isinf(data(:)));
    n = size(sdata, 1);
    if n == 0
      % If a divide by zero error occured, 
      % give it NaN value instead of an empty matrix.
      s = NaN;
    else
      s = mean(sdata, 1);
    end
  else
    for num=1:size(data, dim)
      % Otherwise recurse
      [dims{1:(dim-1)}] = deal(':');
      dims{dim} = num;
      [s(dims{:}) n(dims{:})] = recmean(data(dims{:}), dim - 1);
    end
  end

