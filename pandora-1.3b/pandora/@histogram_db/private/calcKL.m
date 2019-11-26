function [kl_bits] = calcKL(left_data, right_data)

% $Id: calcKL.m 1334 2012-04-19 18:02:13Z cengique $
%
% Author: Cengiz Gunay <cgunay@emory.edu>, 2009/04/03

% Copyright (c) 2007-2009 Cengiz Gunay <cengique@users.sf.net>.
% This work is licensed under the Academic Free License ("AFL")
% v. 3.0. To view a copy of this license, please look at the COPYING
% file distributed with this software or visit
% http://opensource.org/licenses/afl-3.0.php.

% Calculate KL divergence
kl_bits = 0;
for bin_num = 1:size(left_data, 1)
  if left_data(bin_num) ~= 0 && right_data(bin_num) ~= 0
    kl_bits = ...
        kl_bits + ...
        left_data(bin_num) * log2(left_data(bin_num) / ...
                                  right_data(bin_num));
  end
end
