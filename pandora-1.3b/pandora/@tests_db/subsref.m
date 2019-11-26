function b = subsref(a,index)

% subsref - Defines indexing for tests_db objects for () and . operations. 
%
% Usage:
% obj = obj(rows, tests)
% obj = obj.attribute
%
% Description:
% Returns attributes or selects the given test columns and rows
% and returns in a new tests_db object.
%
%   Parameters:
%	obj: A tests_db object.
%	rows: A logical or index vector of rows. If ':', all rows.
%	tests: Cell array of test names or column indices. If ':', all tests.
%	attribute: A tests_db class attribute.
%		
%   Returns:
%	obj: The new tests_db object.
%
% See also: subsref, tests_db
%
% $Id: subsref.m 1335 2012-04-19 18:04:32Z cengique $
%
% Author: Cengiz Gunay <cgunay@emory.edu>, 2004/09/17

% Copyright (c) 2007 Cengiz Gunay <cengique@users.sf.net>.
% This work is licensed under the Academic Free License ("AFL")
% v. 3.0. To view a copy of this license, please look at the COPYING
% file distributed with this software or visit
% http://opensource.org/licenses/afl-3.0.php.

% If a is an array, use built-in methods
if length(a) > 1
  b = builtin('subsref', a, index);
  return;
end

if size(index, 2) > 1
  first = subsref(a, index(1));
  % tail recursion converted to loop
  b = first;
  for index_num = 2:length(index)
    b = subsref(b, index(index_num));
  end
else
  switch index.type
    case '()'
      if length(index.subs) == 1 % Only rows
	b = onlyRowsTests(a, index.subs{1}, ':');
      elseif length(index.subs) == 2 % Both rows and cols
	b = onlyRowsTests(a, index.subs{1:2});
      elseif length(index.subs) == 3 % Rows, cols, and pages.
	b = onlyRowsTests(a, index.subs{1:3});
      else
	error('Only up to three indices are supported: db(rows, cols, pages)');
      end
    case '.'
      b = get(a, index.subs); 
    case '{}'
      b = a{index.subs{:}};
  end
end