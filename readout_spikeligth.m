experimentator= 'SW';
rdata_dir         = 'D:\Postdoc_Margrie\Projects\Callosal\invitro_ephys\Organized';%data directory of raw data;change accordingly
adata_dir         = 'D:\Postdoc_Margrie\Projects\Callosal\output';%data directory of extracted date;change accordingly
%adata_dir         = 'D:\Postdoc_Margrie\Projects\Whisker\output_structure';%data directory of extracted date;change accordingly
ExpXls            = 'C:\Users\simonw\S1-V1 interaction Dropbox\Simon Weiler\Callosal_L6\slice_ephys_structure\Experiment_list.xlsx';%directory where excel batch file is located;change accordingly
%ExpXls            = 'C:\Users\simonw\S1-V1 interaction Dropbox\Simon Weiler\Callosal_L6\slice_ephys_structure\Experiment_list_S1V1.xlsx'
%% parse Experiments XLS database
batchopt          = parseExperimentsXls_ephys(ExpXls);%calls the nested function parseExperimentsXls_ephys and considers the user flag (1 or 0)
nummice           = length(batchopt.mouse);%length of experiments to be analyzed

adder=1;%counting variable
 for i=1:nummice%for loop over experiments across days
     datapath=fullfile(rdata_dir, num2str(batchopt.mouse{i}), filesep);%directory and name of experiments (from excel sheet)
     cd(char(datapath));%go to directory
      for k=1:length(batchopt.exp_ids{i})%loop in bigger loop for each cell per experimental day
          if batchopt.exp_ids{i}(k)<10%for cells with id less then XX0010, e.g., XX0001-XX0009
                n_str = sprintf( '%04d', batchopt.exp_ids{i}(k));
          else
                n_str = sprintf( '%04d', batchopt.exp_ids{i}(k));%for cells with id mor then XX0010, e.g., XX0011-XX0099
          end
          fold_name=[experimentator n_str];%complete cell folder name such as SW0001 or MF0001
            exp_folder=fullfile(datapath,fold_name);%complete folder and directory
            cell_ab=[num2str(batchopt.mouse{i}) fold_name];
          % Check how many files are in the cell folder and cd to it
            list=dir([char(exp_folder) '\*.h5']);
            cd(exp_folder);
            
      end
 end