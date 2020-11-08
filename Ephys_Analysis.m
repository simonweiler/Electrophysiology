%Analysis of slice electropysiology Margrie Lab started Spetember 2019
%% NOTES
%Data information is in excel sheet which can be called by Matlab 
%Nested functions
   % parseExperimentsXls_ephys
   % loadDataFile_wavesurfer
   % Pandora Toolbox 
   % spike event
   
%% Scripts starts here
%% 
LED_stim=1;
intr_prop=1;
image_prop=0;
savefile=1;


%% Set directories and experimentator
experimentator= 'SW';
rdata_dir         = 'D:\Postdoc_Margrie\Projects\Callosal\invitro_ephys\Organized';%data directory of raw data;change accordingly
adata_dir         = 'D:\Postdoc_Margrie\Projects\Callosal\output';%data directory of extracted date;change accordingly
ExpXls            = 'D:\Postdoc_Margrie\Projects\Callosal\experiment_list_slice\Experiment_list.xlsx';%directory where excel batch file is located;change accordingly

%% parse Experiments XLS database
batchopt          = parseExperimentsXls_ephys(ExpXls);%calls the nested function parseExperimentsXls_ephys and considers the user flag (1 or 0)
nummice           = length(batchopt.mouse);%length of experiments to be analyzed
 %% Go through each mouse and recordings iteratively
    
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
            
           %Check wether intrinsic properties or LED stimulation was
           %used 
            idx_led=find(contains({list(:).name},'LED')==1);
            
            intrlist={'IV','Passive','Rheobase','Passive','Sag','Ramp'};
            idx_intr=find(contains({list(:).name},intrlist)==1);
            %check which amplifier
            ampli=batchopt.amplifier{i} (k);
            
              %setting up structure
            Ephys(adder).animal_name=[char(batchopt.mouseID{i})];
            Ephys(adder).patching_date=batchopt.mouse{i};
            Ephys(adder).celID=fold_name;
            Ephys(adder).slice=batchopt.slicen{i} (k);
            Ephys(adder).geno=batchopt.geno{i} (k);
            Ephys(adder).sag=batchopt.sag{i} (k);
            Ephys(adder).label=batchopt.labelcell{i} (k);
            Ephys(adder).area=batchopt.brainarea{i} (k);
            Ephys(adder).wc=batchopt.wholecell{i} (k);
            Ephys(adder).sol=batchopt.solution{i} (k);
            Ephys(adder).amp=batchopt.amplifier{i} (k);
            Ephys(adder).pair=batchopt.paired{i} (k);
            Ephys(adder).distx=batchopt.distancex{i} (k);
            Ephys(adder).disty=batchopt.distancey{i} (k);
            
        if LED_stim==1
            %Get LED stimuli of the cell
            if isempty(idx_led)==1
               disp('No LED evoked stimulus measured');
            else isempty(idx_led)==0
                count1=1;
                count2=1;
                count3=1;
                count4=1;
                for m=1:length(idx_led);  
                filename=[char(exp_folder) '\' list(idx_led(m)).name];
                %load file using loadDataFile_wavesurfer
                data = loadDataFile_wavesurfer(filename);
%                
                stimuli_type=data.header.DataFileBaseName;
                if data.header.AIChannelUnits{1}=='pA';
                    clamp=0;
                   
                else data.header.AIChannelUnits{1}=='mV';
                    clamp=1;
                end
                
                %call function ephys_LED
                if contains(stimuli_type,'train')==1 || contains(stimuli_type,'square')==1   
                [train_n(:,count1) train_p(:,count1) dpeak_n(:,count1) dpeak_p(:,count1) sub_traces_train(:,count1)] = ephys_LED(filename,data,clamp,ampli);
                clamp_v(count1)=clamp;
                count1=count1+1;
                elseif contains(stimuli_type,'ladder')==1;
                    [ladder_n(:,count2) ladder_p(:,count2)  dpeak_n_ladder(:,count2) dpeak_p_ladder(:,count2) sub_traces_ladder(:,count2)] = ephys_LED(filename,data,clamp,ampli);
                clamp_v(count2)=clamp;
                    count2=count2+1;
                elseif contains(stimuli_type,'high')==1;
                    if contains(stimuli_type,'high2')==1;
                   [high_n(:,count3) high_p(:,count3) dpeak_n_high(:,count3) dpeak_p_high(:,count3) sub_traces_high(:,count3)] = ephys_LED(filename,data,clamp,ampli);
                clamp_v(count3)=clamp;
                    count3=count3+1;  
                    else
                     [highf_n(:,count3) highf_p(:,count3) dpeak_n_highf(:,count3) dpeak_p_highf(:,count3) sub_traces_highf(:,count3)] = ephys_LED(filename,data,clamp,ampli);
                      clamp_v(count4)=clamp;
                    count4=count4+1;     
                    end
                else
                    disp('fail');
                end
                
                end
               
            end
            
            Ephys(adder).train_n=train_n;
            Ephys(adder).train_p=train_p;
            Ephys(adder).dpeak_n=dpeak_n;
            Ephys(adder).dpeak_p=dpeak_p;
            Ephys(adder).sub_traces_train=sub_traces_train;
            Ephys(adder).clamp=clamp_v;
           
            if exist('ladder_n')==1;
                Ephys(adder).ladder_n=ladder_n;
                Ephys(adder).ladder_p=ladder_p;
                Ephys(adder).sub_traces_ladder=sub_traces_ladder;
            else
                Ephys(adder).ladder_n=[];
                Ephys(adder).ladder_p=[];
                Ephys(adder).sub_traces_ladder=[];
            end
           
            
            if exist('high_n')==1;
                Ephys(adder).high_n=high_n;
                Ephys(adder).high_p=high_p;
                Ephys(adder).sub_traces_high=sub_traces_high;
            else
                Ephys(adder).high_n=[];
                Ephys(adder).high_p=[];
                Ephys(adder).sub_traces_high=[];
            end
           
              if exist('highf_n')==1;
                Ephys(adder).highf_n=highf_n;
                Ephys(adder).highf_p=highf_p;
                Ephys(adder).sub_traces_highf=sub_traces_highf;
            else
                Ephys(adder).highf_n=[];
                Ephys(adder).highf_p=[];
                Ephys(adder).sub_traces_highf=[];
            end
            
             count1=[];
             count2=[];
             count3=[];
             count4=[];
             train_n=[];
             train_p=[];
             dpeak_n=[];
             dpeak_p=[];
             sub_traces_train=[];
             ladder_n=[];
             ladder_p=[];
             dpeak_n_ladder=[];
             dpeak_p_ladder=[];
             sub_traces_ladder=[];
             high_n=[];
             high_p=[];
             dpeak_n_high=[];
             dpeak_p_high=[];
             sub_traces_high=[];
             highf_n=[];
             highf_p=[];
             dpeak_n_highf=[];
             dpeak_p_highf=[];
             sub_traces_highf=[];
             
             if intr_prop==0;
                 adder=adder+1;
             end
        end
                
                
                
     %% Intrinsic properties
       if intr_prop==1;
            %Get the intrinsic properties of the cell
            if isempty(idx_intr)==1
                disp('No intrinsic properties measured');
                IV=[];
                Rheobase=[];
                Passive=[];
                Sag=[];
                Ramp=[];
                
            Ephys(adder).IV=IV;
            Ephys(adder).Rheobase=Rheobase;
            Ephys(adder).Passive=Passive;
            Ephys(adder).Sag=Sag;
            Ephys(adder).Ramp=Ramp;
            adder=adder+1; 
            else isempty(idx_intr)==0;
                
             [IV Rheobase Passive Sag Ramp] = ephys_intr(list, idx_intr, exp_folder); 
               
            Ephys(adder).IV=IV;
            Ephys(adder).Rheobase=Rheobase;
            Ephys(adder).Passive=Passive;
            Ephys(adder).Sag=Sag;
            Ephys(adder).Ramp=Ramp;
            
            adder=adder+1;      
                    
            end
       end
       %% Image L6 pial depth read out
       if image_prop==1;
       [dis] = img_prop(list,exp_folder);
%        Ephys(adder).pial_depth=pial_depth;
%        Ephys(adder).xi=xi;
%        Ephys(adder).yi=yi;
       Ephys(adder).dis=dis;
       adder=adder+1; 
       end
       
            %% 
             
%             if exist('ramp')==1;
%                 Ephys(adder).ladder_n=ladder_n;
%                 Ephys(adder).ladder_p=ladder_p;
%                 Ephys(adder).sub_traces_ladder=sub_traces_ladder;
%             else
%                 Ephys(adder).ladder_n=[];
%                 Ephys(adder).ladder_p=[];
%                 Ephys(adder).sub_traces_ladder=[];
%             end
           
            
            end
      end
       
      
      
       % SAVE in analyzed directory
    if savefile==1
        cd(adata_dir);
        FileName=['Data_',experimentator,'_',datestr(now, 'hh-dd-mmm-yyyy')];
        %         save(FileName,'-struct','LGN','-v7.3');
        save(FileName,'Ephys','-v7.3');
        disp('FILE SAVED');
    else
        disp('FILE NOT SAVED');
    end
