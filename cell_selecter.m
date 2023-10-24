function filter_out = cell_selecter(data,  varargin)


p = inputParser();
p.KeepUnmatched = false;
p.CaseSensitive = false;
p.StructExpand  = false;
p.PartialMatching = true;

% filters based on data structure
addParameter(p, 'label', nan); 
addParameter(p, 'area' , nan) % when the filter is set to nan it is not used 
addParameter(p, 'wc', nan);
addParameter(p, 'sol', nan); 
addParameter(p, 'amp', nan);
addParameter(p, 'pair', nan);
addParameter(p, 'distx', nan);
addParameter(p, 'disty', nan);
addParameter(p, 'geno', nan);
addParameter(p, 'drugs', nan);
addParameter(p, 'optovariant', nan);
addParameter(p, 'layer', nan);
addParameter(p, 'qualityinput', nan);
% % filters based on QC sheet
% addParameter(p, 'transduction_QC', 1); % 
% addParameter(p, 'Conf_QC', nan);
% addParameter(p, 'morph_QC', nan);
% addParameter(p, 'TwoPTrans_QC', nan);
% addParameter(p, 'CCF_QC', nan);
% addParameter(p, 'A_ramp_signal_QC', nan);
% addParameter(p, 'N_ramp_signal_QC', nan);
% addParameter(p, 'AtoN_Rschange_QC', nan);
% addParameter(p, 'CS_internal', nan);
% addParameter(p, 'Interneuron_morph', 0); % should be 0

parse(p, varargin{:})
cell_filters = p.Results;

% remove fields with nan
fn = fieldnames(cell_filters);
for k=1:numel(fn)
    if isnumeric(cell_filters.(fn{k})) && isnan(cell_filters.(fn{k}))
        cell_filters = rmfield(cell_filters,fn{k});
    end
end

% % get QC sheet
% basedir = 'C:\Users\slice setup\SimonW\LGN\';
% originaldir = cd;
% cd([basedir])
% d = dir('QC_Filter*.xlsx');
% file_name = d(find([d(:).datenum]==max([d(:).datenum]))).name;
% QC_table = readtable([file_name]);
% QC_table = table2struct(QC_table);
% cd(originaldir) % go back to where you came from!
% 
% data_cell_names = join([cellfun(@(x) x(1:6),{data(:).patching_date},'UniformOutput',false)' ...
%     {data(:).experimentator}'...
%     {data(:).cellname}']);
% QC_table_cell_names = join([{QC_table(:).file_name}' ...
%     {QC_table(:).experimentator}'...
%     {QC_table(:).cellname}']);

% % remove cells not present in data from QC_table
% try
%     QC_table(~ismember(unique(flip(QC_table_cell_names),unique(data_cell_names)))) = [];
% end
% 
% % match order of QC_table and add QC subfields to data (do not overwrite those already present!)
% I = cell2mat(cellfun(@(x) strmatch(x,data_cell_names),QC_table_cell_names, 'UniformOutput', false));
% fn_QC = fieldnames(QC_table);
% fn_data = fieldnames(data);
% data2 = data; 
% for fn_idx = 6:length(fn_QC)
%     for celln_QC = 1:length(I)
%         celln_data = I(celln_QC);
%         data2(celln_data).(fn_QC{fn_idx}) = QC_table(celln_QC).(fn_QC{fn_idx});
%     end
%     for celln_data = 1:length(data)
%         if isempty(data2(celln_data).(fn_QC{fn_idx}))
%             data2(celln_data).(fn_QC{fn_idx}) = nan;
%         end
%     end
% end
data2 = data;
% loop through filter and create filter vectors
fn = fieldnames(cell_filters);
filters_array = [];
for k = 1:length(fn)
    if strcmp(fn{k},'amp')
        temp_filter = cellfun(@(x) all(cell_filters.amp==x),{data2.(fn{k})});
        filters_array(k,:) = temp_filter;
    elseif strcmp(fn{k},'distx')
        temp_filter = cell_filters.transduction_QC <= [data2.transduction_QC];
        filters_array(k,:) = temp_filter;
    elseif strcmp(fn{k},'disty')
        temp_filter = cellfun(@(x) ~isempty(x),{data2.scracm_red});
        filters_array(k,:) = temp_filter;
    elseif strcmp(fn{k},'label')
        temp_filter = cellfun(@(x) all(cell_filters.label==x),{data2.(fn{k})});
        filters_array(k,:) = temp_filter;
   
     elseif strcmp(fn{k},'sol')
       temp_filter = cellfun(@(x) all(cell_filters.sol==x),{data2.(fn{k})});
        filters_array(k,:) = temp_filter;
         elseif strcmp(fn{k},'geno')
       temp_filter = cellfun(@(x) all(cell_filters.geno==x),{data2.(fn{k})});
        filters_array(k,:) = temp_filter;
           elseif strcmp(fn{k},'area')
       temp_filter = cellfun(@(x) all(cell_filters.area==x),{data2.(fn{k})});
        filters_array(k,:) = temp_filter;
        elseif strcmp(fn{k},'wc')
       temp_filter = cellfun(@(x) all(cell_filters.wc==x),{data2.(fn{k})});
        filters_array(k,:) = temp_filter;
    
          elseif strcmp(fn{k},'pair')
       temp_filter = cellfun(@(x) all(cell_filters.pair==x),{data2.(fn{k})});
        filters_array(k,:) = temp_filter;
              elseif strcmp(fn{k},'drugs')
       temp_filter = cellfun(@(x) all(cell_filters.drugs==x),{data2.(fn{k})});
        filters_array(k,:) = temp_filter;
            elseif strcmp(fn{k},'optovariant')
       temp_filter = cellfun(@(x) all(cell_filters.optovariant==x),{data2.(fn{k})});
        filters_array(k,:) = temp_filter;
          elseif strcmp(fn{k},'layer')
       temp_filter = cellfun(@(x) all(cell_filters.layer==x),{data2.(fn{k})});
        filters_array(k,:) = temp_filter;
           elseif strcmp(fn{k},'qualityinput')
       temp_filter = cellfun(@(x) all(cell_filters.qualityinput==x),{data2.(fn{k})});
        filters_array(k,:) = temp_filter;
    else
        temp_filter = cellfun(@(x) all(cell_filters.(fn{k})==x),{data2.(fn{k})});
        filters_array(k,:) = temp_filter;
    end
end
if isempty(filters_array)
    filter_out = ones(1,length(data2));
else
    filter_out = all(filters_array,1);
end
disp([num2str(sum(filter_out)) ' cells pass all filters' ])
    
    