function [data, time_trace] = readgenbin(filename, start_time, end_time)
%
% READGENESIS -- Read data from a binary GENESIS file
%
% Usage:
%
%   [data, time_trace] = readgenbin(filename, start_time, end_time);
%
%   filename             -- no checking for binary type is made, so if you
%                           want reliability please ensure the file is a binary.
%   
%   start_time, end_time -- time in milliseconds relative to the ACTUAL 
%                           time of the experiment at which data adquisition
%                           started (if you start gathering data at 200 ms
%                           and you specify 0 start time it will not work).
%
%   data                  -- data return.
%   time_trace            -- corresponding time (in ms) for data.
%
% <adelgado@biology.emory.edu> based in open binary from Simon Peron.
%

% Open the file
fil = fopen(filename, 'r');

% Figure out how many bites we can read at most
fseek(fil, 0, 'eof');
length = ftell(fil);

% Determine the sampling freq. and number of channels from header
fseek(fil, 80, 'bof');
act_start_time = fread(fil, [1 1], 'float32');
timestep = fread(fil, [1 1], 'float32');
freq = 1/timestep;
adjustment_ratio = freq * 1e-3; 
num_cols = fread(fil, [1 1], 'int');
ev_start_time = round(start_time * adjustment_ratio);
ev_end_time = round(end_time * adjustment_ratio);

% Check to make sure the user did not specify invalid time ranges
header_size = 96 + (3 * num_cols * 4);
read_length = (ev_end_time - ev_start_time);
read_start = header_size + (ev_start_time * 4 * num_cols);

% Read those data in
fseek(fil, round(read_start), 'bof');
A = fread(fil, [num_cols read_length], 'float32');
fseek(fil, round(read_start), 'bof');
fclose(fil);

data = A';

% Create the time_trace
if nargout > 1
    time_size = end_time - start_time;
    for i = 1:(round(time_size*adjustment_ratio))
	time_trace(i, 1) = (i/adjustment_ratio)+start_time;
    end
end

% clear ev_*, act_*, num_ % they might want them around
