function obj = trace(data_src, dt, dy, id, props)

% trace - Load a generic data trace. It can be membrane voltage, current, etc.
%
% Usage:
% obj = trace(data_src, dt, dy, id, props)
%
% Description:
%   This object is designed to recognize most data file formats. See the
% data_src parameter below. Traces for specific experimental or simulation
% protocols can extend this class for adding new parameters.
%
% Parameters:
%    data_src: Trace data as a column vector OR name of a data file generated by either 
%	Genesis (.bin, .gbin, .genflac), Neuron, PCDX (.all), Matlab (.mat) or
%	NeuroShare (.nsn, .nev, .stb, .plx, .nex, .map, .son, .smr, .mcd) files.
%    dt: Time resolution in [s], unless specified in HDF5 or NeuroShare
%    	file. For Neuron ASCII files, it is used as the file data units.
%    dy: y-axis resolution in [ISI (V, A, etc.)], unless specified in HDF5 or NeuroShare file.
%    id: Identification string
%    props: A structure with any optional properties.
%	  scale_y: Y-axis scale to be applied to loaded data.
%	  offset_y: Y-axis offset to be added to loaded and scaled data.
%	  unit_y: Unit of Y-axis as in 'V' or 'A' (default='V').
%         y_label: String to put on Y-axis of plots.
%	  trace_time_start: Samples in the beginning to discard [dt]
%	  baseline: Resting potential.
%	  channel: Channel to read from file Genesis, PCDX, NeuroShare or
%	  	   Neuron file, or column in a data vector.
%         numTraces: Divide the single column vector of data into this
%         		many columns by making it a matrix.
%	  file_type: Specify file type instead of guessing from extension:
%		'genesis': Raw binary files created with Genesis disk_out method.
%		'genesis_flac': Compressed Genesis binary files.
%		'neuron': Binary files created with Neuron's Vector.vwrite method.
%		'neuronascii': Ascii files created from Neuron's Vector objects. 
%			       Uses time step in file to scale given dt (Must be in ms).
%		'pcdx': .ALL data acquisition files from PCDX program.
%		'matlab': Matlab .MAT binary files with matrix data.
%		'neuroshare': One of the vendor formats recognized by
%			NeuroShare Windows DLLs. See above and http://neuroshare.org. A
%			scale_y value may need to be supplied to get the correct units.
%		'abf': AxoClamp .ABF format read with abf2load from Matlab FileExchange.
%         file_endian: 'l' for little endian and 'b' for big endian
%           (default='n', for native endian). See machineformat option 
%	    in fopen for more info.
%         traces: Trace numbers as a numeric array or as a string with
%         	numeric ranges (e.g., '1 2 5-10 28') for PCDX files.
%	  spike_finder: Method of finding spikes 
%	                (1 for findFilteredSpikes, 2 for Li Su's
%	                findspikes, and 3 for Alfonso Delgado Reyes's 
%			findspikes_old). Methods 2 and 3 require a threshold.
%	  threshold: Spike finding threshold. For the findspikes method,
%	  	     it is either a scalar, or [thres1 thres2] to define
%	  	     a range. For findFilteredSpikes it is used on the
%	  	     filtered data and the default is 2/3 max amplitude
%	  	     of band-passed data, but with a minimum of 15. 
%         downThreshold: (Only for findFilteredSpikes) Size of the trough
%         	     after the spike peak in filtered data (Default=-2).
%	  minInit2MaxAmp, minMin2MaxAmp: For spike_shape elimination,
%	  	     conditions of minimal allowed values for 
%	  	     initial point to max point and minimal point to 
%	             max point, respectively (Default=10 for both).
%	  init_Vm_method: Method of finding spike thresholds during spike
%	  		shape calculation (see spike_shape/spike_shape).
%	  init_threshold: Spike initiation threshold (deriv or accel).
%			(see above methods and implementation in calcInitVm)
%	  init_lo_thr, init_hi_thr: Low and high thresholds for slope.
%         custom_filter: Recommended if sampling rate differs appreciably from 10 kHz.
%                      If custom_filter == 1, a filter with custom lowpass and highpass
%                      cutoffs can be specified. This allows for fast and accurate spike
%                      discrimination. The filter type used is a 2-pole
%                      butterworth, different than the default high-order
%                      Cheby2. Creates new prop called 'butterWorth' to
%                      hold the filter.
%         lowPassFreq: If set, it sets a new low pass cutoff for custom filter. Default is 3000Hz
%         highPassFreq: If set it sets a new high pass cutoff for custom filter. Default is 50 Hz
%	  quiet: If 1, reduces the amount of textual description in plots
%	  	and does not add information to id field.
%		
%   Returns a structure object with the following fields:
%	data: The trace column matrix.
%	dt, dy, id, props (see above)
%
% General operations on trace objects:
%   trace		- Construct a new trace object.
%   plot		- Graph the trace.
%   display		- Returns and displays the identification string.
%   get			- Gets attributes of this object and parents.
%   subsref		- Allows usage of . operator.
%   calcAvg		- Calculate the average value of the trace period.
%   calcMin		- Calculate the minimum value of the trace period.
%   calcMax		- Calculate the minimum value of the trace period.
%   periodWhole		- Return the whole period of this trace.
%   findFilteredSpikes	- Use a band-pass filter to smooth the data and
%			find spikes afterwards. 
%   getResults		- Calculates a set of tests.
%   spike_shape		- Build a trace of the average spike shape in here.
%   spikes		- Build a spikes object with the spikes found here.
%
% Converter methods:
%   spikes		- Find the spikes and construct a spikes object.
%   spike_shape		- Construct a spike_shape object from this trace.
% 
% Additional methods:
%	See methods('trace')
%
% See also: spikes, spike_shape, cip_trace, period
%
% $Id: trace.m 1354 2013-02-06 19:16:02Z cengique $
%
% Author: Cengiz Gunay <cgunay@emory.edu>, 2004/07/30
% Modified:
%   - allow custom filter, Thomas D. Sangrey 2007/12/04

% Copyright (c) 2007 Cengiz Gunay <cengique@users.sf.net>.
% This work is licensed under the Academic Free License ("AFL")
% v. 3.0. To view a copy of this license, please look at the COPYING
% file distributed with this software or visit
% http://opensource.org/licenses/afl-3.0.php.

vs = warning('query', 'verbose');
verbose = strcmp(vs.state, 'on');

if nargin == 0 % Called with no params
   obj.data = [];
   obj.dt = 1;
   obj.dy = 1;
   obj.id = '';
   obj.props = struct([]);
   obj = class(obj,'trace');
 elseif isa(data_src,'trace') % copy constructor?
   obj = data_src;
else
   if ~ exist('props','var')
     props = struct;
   end

   if ~ exist('id','var')
     id = '';
   end

   if isa(data_src, 'char') % filename?
     [path, filename, ext] = fileparts(data_src);

     ext = lower(ext); % Case insensitive matches for file extension

     % Default if native endian of this computer
     if ~ isfield(props, 'file_endian')
       [c,maxsize,endian] = computer;
       props.file_endian = endian;
     end
     
     % if file type not specified, use file extension to guess it
     if ~ isfield(props, 'file_type')
       if strcmpi(ext, '.bin') || strcmpi(ext, '.gbin') % Genesis file
         props.file_type = 'genesis';
       else
         props.file_type = '';
       end
     end

     if strcmpi(props.file_type, 'genesis') 
       channel = 1; % by default
       if isfield(props, 'channel')
        channel = props.channel;
       end

       if ~ isempty(findstr(filename, '_BE_')) || ...
           ~ isempty(findstr(filename, '_BE.')) || ...
           strcmpi(props.file_endian, 'b')
         % check if readgenesis.mex* file available
         if exist('readgenesis_BE') == 3
           % Use big-endian (Mac, Sun) version of readgenesis
           data = readgenesis_BE(data_src, channel);
         else
           % otherwise use native Matlab reader
           data = readgenbin(data_src, NaN, NaN, props.file_endian);
           data = data(:, channel);
         end
       else
         % check if readgenesis.mex* file available
         if exist('readgenesis') == 3
           % Use regular (i386 PCs) little-endian version of readgenesis
           data = readgenesis(data_src, channel);
         else
           % otherwise use native Matlab reader
           data = readgenbin(data_src, NaN, NaN, props.file_endian);
           data = data(:, channel);
         end           
       end

     elseif strcmpi(props.file_type, 'genesis_flac') || ...
	   strcmpi(ext, '.genflac') % Compressed 16-bit genesis file
       channel = 1; % by default
       if isfield(props, 'channel')
         channel = props.channel;
       end
       data = readgenesis16bit(data_src);
       data = data(:, channel);

     elseif strcmpi(props.file_type, 'neuron')
       [data err] = readNeuronVecBin(data_src, props.file_endian);
       if isempty(strfind(err, 'Success'))
         error([ 'Failed to load Neuron binary ''' data_src ''' because of ' ...
                 'error: ' err ]);
       end
       channel = ':'; % by default
       if isfield(props, 'channel')
         channel = props.channel;
       end
       data = data(:, channel);
     
     elseif strcmpi(props.file_type, 'neuronascii')
       [data, label] = readNeuronVecAscii(data_src);
       dt = dt * diff(data(1:2, 1)); % scale given dt by delta in first column
       data = data(:, 2);
       id = [ id label];

     elseif strcmpi(props.file_type, 'pcdx') || ...
	   strcmpi(ext, '.all') % PCDX file
       %disp('Loading PCDX trace');
       data = loadtraces(data_src, props.traces, props.channel, 1);
       
     elseif strcmpi(props.file_type, 'matlab') || ...
	   strcmpi(ext, '.mat') % MatLab file
       s = load(data_src);
       fields = fieldnames(s);
       data = getfield(s, fields{1});	% Assuming there's only one vector
       
     elseif strcmpi(props.file_type, 'hdf5') || ...
	   strcmpi(ext, '.hdf5') 
       % new neurosage file
       % reset file type explicitly
       props.file_type = 'hdf5';
       s1 = ns_read(props.AcquisitionData{props.channel});
        % Make sure the data is in a column vector
        if size(s1.Y, 1) > size(s1.Y, 2)
          data = s1.Y;
        else
          data = s1.Y';
        end
        % NeuroShare files with these extensions:
        % nsn, .nev, .stb, .plx, .nex, .map, .son, .smr, .mcd
     elseif strcmpi(props.file_type, 'neuroshare') || ...
         strcmpi(ext, '.nsn') || strcmpi(ext, '.nev') || strcmpi(ext, '.stb') || ...
         strcmpi(ext, '.plx') || strcmpi(ext, '.nex') || strcmpi(ext, '.map') || ...
         strcmpi(ext, '.son') || strcmpi(ext, '.smr') || strcmpi(ext, ...
                                                         '.mcd')
       [data, dt, id, props] = trace_loadns(data_src, props);
       
       % assume units in mV? Nothing specified in Matlab loader
       if isempty(dy), dy = 1e-3; end
     elseif strcmpi(ext, '.abf') || strcmpi(props.file_type, 'abf')
       % Load AxoClamp ABF format using MathWorks fileExchange entry
       % abf2load
       [dt, data, y_units, dy, cell_name] = loadABF(data_src, props);

       if isfield(props, 'channel')
         data = data(:, props.channel);
       end

       props.unit_y = y_units;
       if ~ isfield(props, 'quiet')
         id = [ cell_name id ];
       end

     else
       warning(['No matching load function found for file ''' data_src ...
              ''' or specified type ''' ...
	      props.file_type '''; defaulting to Matlab load() function.']);
       s = load(data_src);
       data = getfield(s, fields{1});	% Assuming there's only one vector
       data = data(:, getFieldDefault(props, 'channel', ':'));
     end

     % use the filename as id unless otherwise specified
     if ~ exist('id','var') || isempty(id)
       id = name;
     end

   elseif isnumeric(data_src)
     data = data_src;
     if size(data_src, 2) > size(data_src, 2)
       data = data';            % convert to column vector
     end
     channel = getFieldDefault(props, 'channel', ':');
     if ~ strcmp(channel, ':')
       data = data(:, channel);
       id = [id ', chan ' num2str(channel) ];
     end
   else
     error('Unrecognized data source!');
   end

   % reshape data if requested
   if isfield(props, 'numTraces')
     total_len = size(data, 1);
     if mod(total_len, props.numTraces) ~= 0
       error(['props.numTraces=' num2str(props.numTraces) ...
              ' does not evenly divide data vector of size=' ...
              num2str(size(data, 1))]);
     end
     data = reshape(data, total_len/props.numTraces, props.numTraces);
   end
   
   % Scale the loaded data if desired
   if isfield(props, 'scale_y')
     % scale HDF5 data only if gain is unspecified in file metadata (by Li Su)
     if ~( (isfield(props, 'file_type') && strcmpi(props.file_type, 'hdf5')) ) || ...
         (isfield(props, 'AcquisitionData') && ...
          strcmpi( ... 
            props.AcquisitionData{props.channel}.PhysicalDevice, ...
            'Unspecified' ...
            ) ...
          )
         data = props.scale_y * data;
         if verbose
           warning(sprintf(['Device unspecified. Trace multiplied by given ' ...
                            'gain %g'], props.scale_y ));
         end
      end
   end

   % Apply offset to data if desired (after scaling?)
   if isfield(props, 'offset_y')
     data = data + props.offset_y;
   end

   % Crop the data if desired
   if isfield(props, 'trace_time_start')
%     data =  data(:, props.trace_time_start:end);
     data =  data(props.trace_time_start:end, :);
   end

   % Custom filter props
   if isfield(props, 'custom_filter') && props.custom_filter == 1
     if isfield(props, 'highPassFreq')
       hpf = props.highPassFreq;
     else
       hpf = 50;
     end
     if isfield(props, 'lowPassFreq')
       lpf = props.lowPassFreq;
     else
       lpf =3000;
     end
     % make a new filter based on the given dt
     % multiply by 2 to find nyquist freq
     assert(all(2*[lpf hpf]*dt <= 1), ...
            [ 'For custom filter, highPassFreq and lowPassFreq must be <= ' ...
              'Nyquist freq (' num2str(1/dt/2) ').' ]);
     [b,a] = butter(2,2*hpf*dt,'high');
     butterWorth.highPass = struct('b', b, 'a', a);
     [b,a] = butter(2,2*lpf*dt,'low');
     butterWorth.lowPass = struct('b', b, 'a', a);
     props.butterWorth = butterWorth;
   end
   obj.data = data;
   obj.dt = dt;
   obj.dy = dy;
   obj.id = id;
   obj.props = props;
   obj = class(obj, 'trace');
end

