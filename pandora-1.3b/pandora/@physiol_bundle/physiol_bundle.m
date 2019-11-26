function a_bundle = physiol_bundle(a_cell, props)

% physiol_bundle - The physiology dataset and the DB created from it bundled together.
%
% Usage:
% a_bundle = physiol_bundle(a_cell, props)
%
%   Parameters:
%	a_cell: A cell array that contains the following elements:
%	a_dataset: A cell-enclosed physiol_cip_traceset_fileset object.
%	a_db: The raw params_tests_db object created from the dataset. 
%		It only needs to have the pAcip, pAbias, TracesetIndex, and ItemIndex columns.
%	a_joined_db: The one-treatment-per-line DB created from the raw DB.
%	props: A structure with any optional properties.
%	  controlDB: Use this as the ontrol DB rather than computing.
%		
% Description:
%   This is a subclass of dataset_db_bundle, specialized for physiology datasets. 
%
% Returns a structure object with the following fields:
%	dataset_db_bundle, 
%	joined_control_db: DB of control neurons (no pharmacological applications).
%
% General operations on physiol_bundle objects:
%   physiol_bundle 	- Construct a new physiol_bundle object.
%   display		- Returns and displays the identification string.
%   get			- Gets attributes of this object and parents.
%   subsref		- Allows usage of . operator.
%
% Additional methods:
%	See methods('physiol_bundle')
%
% See also: dataset_db_bundle, tests_db, params_tests_dataset
%
% $Id: physiol_bundle.m 1335 2012-04-19 18:04:32Z cengique $
%
% Author: Cengiz Gunay <cgunay@emory.edu>, 2005/12/13

% Copyright (c) 2007 Cengiz Gunay <cengique@users.sf.net>.
% This work is licensed under the Academic Free License ("AFL")
% v. 3.0. To view a copy of this license, please look at the COPYING
% file distributed with this software or visit
% http://opensource.org/licenses/afl-3.0.php.

if nargin == 0 % Called with no params
  a_bundle.joined_control_db = params_tests_db;
  a_bundle = class(a_bundle, 'physiol_bundle', dataset_db_bundle);
elseif isa(a_cell, 'physiol_bundle') % copy constructor?
  a_bundle = a_cell;
else
  if ~ exist('props', 'var')
    props = struct([]);
  end

  [a_dataset, a_db, a_joined_db] = deal(a_cell{:});
  
  if ~ isfield(props, 'controlDB')
    a_bundle.joined_control_db = ...
        a_joined_db(a_joined_db(:, {'TTX', 'Apamin', 'EBIO', ...
                        'XE991', 'Cadmium', 'drug_4AP'}) == zeros(1, 6), ...
                    :);
  else
    a_bundle.joined_control_db = props.controlDB;
  end
  
  % TODO: do not remove anything from a_db, it's already stripped
  % used to keep only: (:, {'pAcip', 'pAbias', 'TracesetIndex', 'ItemIndex'})
  a_bundle = ...
      class(a_bundle, 'physiol_bundle', ...
	    dataset_db_bundle(a_dataset, a_db, ...
			      a_joined_db, props));
end

