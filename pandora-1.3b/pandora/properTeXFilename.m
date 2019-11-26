function filename = properTeXFilename(filename)

%  properTeXFilename - Replaces characters in string to make it a valid filename for inclusion in TeX documents.
%
% Usage:
% filename = properTeXFilename( filename )
%
% Parameters:
%	filename: An input filename string (without extension!).
%
% Returns:
% 	filename: The corrected proper filename.
%
% Description:
%   It will replace characters like space, '/', '.', etc.
%
% Example: 
% >> fname = properTeXFilename('hello world/1')
% ans = 'hello_world+1' 
%
% $Id: properTeXFilename.m 1334 2012-04-19 18:02:13Z cengique $
% Author: Cengiz Gunay <cgunay@emory.edu>, 2005/12/20

% Copyright (c) 2007 Cengiz Gunay <cengique@users.sf.net>.
% This work is licensed under the Academic Free License ("AFL")
% v. 3.0. To view a copy of this license, please look at the COPYING
% file distributed with this software or visit
% http://opensource.org/licenses/afl-3.0.php.

%filename = strrep(filename, '(', '_');
%filename = strrep(filename, ')', '_');
%filename = strrep(filename, ';', '_');
filename = strrep(filename, ' ', '_');
filename = strrep(filename, '\_', '_');
filename = strrep(filename, '.', '_');
filename = strrep(filename, '/', '+');  % And the /'s with +'s
