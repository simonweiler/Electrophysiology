function [params, param_names, tests, test_names, docs] = readDBItems(obj, items)

% readDBItems - Reads all items to generate a params_tests_db object.
%
% Usage:
% [params, param_names, tests, test_names] = readDBItems(obj, items)
%
% Description:
%   This is a generic method to convert from params_tests_fileset to
% a params_tests_db, or a subclass. This method depends on the  
% paramNames, testNames, and itemResultsRow functions. 
% Outputs of this function can be directly fed to the constructor of
% a params_tests_db or a subclass.
%
%   Parameters:
%	obj: A params_tests_fileset object.
%	items: (Optional) List of item indices to use to create the db.
%		
%   Returns:
%	params, param_names, tests, test_names: See params_tests_db.
%
% See also: params_tests_db, params_tests_fileset, itemResultsRow
%	    testNames, paramNames
%
% $Id: readDBItems.m 1335 2012-04-19 18:04:32Z cengique $
%
% Author: Cengiz Gunay <cgunay@emory.edu>, 2004/11/24

% Copyright (c) 2007 Cengiz Gunay <cengique@users.sf.net>.
% This work is licensed under the Academic Free License ("AFL")
% v. 3.0. To view a copy of this license, please look at the COPYING
% file distributed with this software or visit
% http://opensource.org/licenses/afl-3.0.php.

if ~ exist('items', 'var')
  items = 1:length(obj.list);
end

try 
  % Collect info for generating the DB
  [test_names a_prof] = testNames(obj, items(1));
  if isa(a_prof, 'params_results_profile')
    param_names = fieldnames(a_prof.params)';
  else
    param_names = paramNames(obj, items(1));
  end
  num_items = length(items);
catch
  err = lasterror;
  warning(['Error caught during database creation, before starting to ' ...
	   'read items: ' err.message '. Returning empty database.']);
  if isfield(err, 'stack')
    disp('Stack trace:');
    for stack_item = 1:length(err.stack)
      disp([ '  ' err.stack(stack_item).file ' at ' num2str(err.stack(stack_item).line) ...
	    ' (' err.stack(stack_item).name ').' ]);
    end
  end
  test_names = [];
  param_names = [];
  params = [];
  tests = [];
  return
end

% Preallocating matrices dramatically speeds up the filling process
params = repmat(0, num_items, length(param_names));
tests = repmat(0, num_items, length(test_names));
docs = cell(num_items, 1);

% Batch process all items
start_time = cputime;

disp('Reading: ');
line_buffer = '';

try 
  row_index = 1;
  doc_row_index = 1;
  for item_index=items
    %disp(sprintf('File number: %d\r', item_index));
    
    % doesn't work on the cluster
    if false && usejava('jvm')
      if checkError(java.lang.System.out)
	disp('error in System.out stream!');
      else
	flush(java.lang.System.out);
	print(java.lang.System.out, [ num2str(item_index) ', ' ]);	
      end
      if mod(row_index, 20) == 0
	disp(' ');
      end
    else
      line_buffer = [line_buffer num2str(item_index) ', ' ];
      if length(line_buffer) > 70
	disp(line_buffer);
	line_buffer = '';
      end
    end
    [params_row, tests_row, a_doc] = itemResultsRow(obj, item_index);

    % allow multiple rows returned by one item
    num_rows = size(tests_row, 1);
    if num_rows > 1
      params(row_index:(row_index+num_rows-1), :) = params_row;
      tests(row_index:(row_index+num_rows-1), :) = tests_row;
      docs{doc_row_index} = a_doc;
      row_index = row_index + num_rows;      
      % TODO: pre-allocate based on new estimate! 
      % e.g., params(est_row, :) = repmat(NaN, ...)
    else
      params(row_index, :) = params_row;
      tests(row_index, :) = tests_row;
      docs{doc_row_index} = a_doc;
      row_index = row_index + 1;
    end
    doc_row_index = doc_row_index + 1;
  end
catch
  err = lasterror;
  % Matlab R2006b has debug file information in err
  warning(['Error caught during database creation at item index ' ...
	   num2str(item_index) ': ' err.message '. Truncating database.']);
  disp([ 'Item in question is:']);
  disp(getItem(obj, item_index));
  if isfield(err, 'stack')
    disp('Stack trace:');
    for stack_item = 1:length(err.stack)
      disp([ '  ' err.stack(stack_item).file ' at ' num2str(err.stack(stack_item).line) ...
	    ' (' err.stack(stack_item).name ').' ]);
    end
  end
  params(row_index:size(params, 1), :) = [];
  tests(row_index:size(tests, 1), :) = [];
end

% put the docs together
docs = doc_multi(docs, get(obj, 'id'));

end_time = cputime;

disp(sprintf('Elapsed time took %.2f seconds.', end_time - start_time));

