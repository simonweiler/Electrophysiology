str_L6    = 'D:\Postdoc_Margrie\Projects\Callosal\output';
folder_list = uipickfiles('FilterSpec',str_L6);
load(char(folder_list));
%sampling rate
srF=20;
sr=20000;
%% 
ntsr_ac=cell_selecter(Ephys,'label',4,'geno',6,'sol',2,'optovariant',1,'area',3);
cp_ac=cell_selecter(Ephys,'label',1,'geno',6,'sol',2,'optovariant',1,'area',3);
%% 
[ntsr_psc] = readout_amp_update(Ephys,ntsr_ac ,1,2,1,2);
[ntsr_psc_cshf] = readout_amp_update(Ephys,ntsr_ac ,2,2,1,2);
[ntsr_psc_cshf2] = readout_amp_update(Ephys,ntsr_ac ,3,2,1,2);
[cp_psc] = readout_amp_update(Ephys,cp_ac ,1,2,1,2);
[cp_psc_cshf] = readout_amp_update(Ephys,cp_ac ,2,2,1,2);
[cp_psc_cshf2] = readout_amp_update(Ephys,cp_ac ,3,2,1,2);
%% 

tempd=[];tempd2=[];
col=3;
tempd=[ntsr_psc(:,col)];
tempd2=[cp_psc(:,col)];
data=[];data=[tempd2 tempd];
paired_plot_box(data,{'g','m'});
xticklabels({'CP','Ntsr1'});
ylabel('Light evoked EPSC amplitude');