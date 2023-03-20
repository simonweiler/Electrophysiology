function [ntsr_cs_cells ntsr_cpncs_cells ntsr_k_cells ntsr_cpnk_cells ...
    penk_cs_cells penk_cpncs_cells penk_k_cells penk_cpnk_cells ...
    cci_k_cells cci_cpnk_cells som_cs_cells som_cpncs_cells som_k_cells som_cpnk_cells ...
    pv_cs_cells pv_cpncs_cells pv_k_cells pv_cpnk_cells vip_k_cells vip_cpnk_cells] = readout_input_strength_callosal(Ephys)
%% compare CPN to NTSR1, SOM and PV EPSC input PAIRED
%CS NTSR1 paired (no drugs and before washin)
temp1=[];ntsr_cs1=[];
for i=1:6
temp1(i,:)=cell_selecter(Ephys,'drugs',0,'label',4,'geno',6,'sol',2,'optovariant',1,'pair',i);
end
ntsr_cs1=sum(temp1);
temp1=[];ntsr_cs2=[];
for i=1:6
temp1(i,:)=cell_selecter(Ephys,'drugs',1,'label',4,'geno',6,'sol',2,'optovariant',1,'pair',i);
end
ntsr_cs2=sum(temp1);
ntsr_cs=[];
ntsr_cs=ntsr_cs1+ntsr_cs2;

%CS NTSR1 paired (no drugs and before washin) CORRESPONDING CPN cells
temp1=[];ntsr_cpncs1=[];
for i=1:6
temp1(i,:)=cell_selecter(Ephys,'drugs',0,'label',1,'geno',6,'sol',2,'optovariant',1,'pair',i);
end
ntsr_cpncs1=sum(temp1);
%
temp1=[];ntsr_cpncs2=[];
for i=1:6
temp1(i,:)=cell_selecter(Ephys,'drugs',1,'label',1,'geno',6,'sol',2,'optovariant',1,'pair',i);
end
ntsr_cpncs2=sum(temp1);
ntsr_cpncs=[];
ntsr_cpncs=ntsr_cpncs1+ntsr_cpncs2;
%remove cell 253 from ntsr1 cause its pair is with TTX
ntsr_cs(253)=0;


%same for K
%CS NTSR1 paired (no drugs and before washin)
temp1=[];ntsr_k1=[];
for i=1:6
temp1(i,:)=cell_selecter(Ephys,'drugs',0,'label',4,'geno',6,'sol',1,'optovariant',1,'pair',i);
end
ntsr_k1=sum(temp1);

temp1=[];ntsr_k2=[];
for i=1:6
temp1(i,:)=cell_selecter(Ephys,'drugs',1,'label',4,'geno',6,'sol',1,'optovariant',1,'pair',i);
end
ntsr_k2=sum(temp1);
ntsr_k=[];
ntsr_k=ntsr_k1+ntsr_k2;

%CS NTSR1 paired (no drugs and before washin) CORRESPONDING CPN cells
temp1=[];ntsr_cpnk1=[];
for i=1:6
temp1(i,:)=cell_selecter(Ephys,'drugs',0,'label',1,'geno',6,'sol',1,'optovariant',1,'pair',i);
end
ntsr_cpnk1=sum(temp1);
%
temp1=[];ntsr_cpnk2=[];
for i=1:6
temp1(i,:)=cell_selecter(Ephys,'drugs',1,'label',1,'geno',6,'sol',1,'optovariant',1,'pair',i);
end
ntsr_cpnk2=sum(temp1);
ntsr_cpnk=[];
ntsr_cpnk=ntsr_cpnk1+ntsr_cpnk2;
%remove cell 259 and 261 from ntsr1 cause its pair is with TTX
ntsr_k([259 261 329])=0;
%% PENK
%CS
temp1=[];penk_cs=[];
for i=1:6
temp1(i,:)=cell_selecter(Ephys,'drugs',0,'label',5,'geno',7,'sol',2,'optovariant',1,'pair',i);
end
penk_cs=sum(temp1);
%K
temp1=[];penk_k=[];
for i=1:6
temp1(i,:)=cell_selecter(Ephys,'drugs',0,'label',5,'geno',7,'sol',1,'optovariant',1,'pair',i);
end
penk_k=sum(temp1);

%Penk CPN
%CS
temp1=[];penk_cpncs=[];
for i=1:6
temp1(i,:)=cell_selecter(Ephys,'drugs',0,'label',1,'geno',7,'sol',2,'optovariant',1,'pair',i);
end
penk_cpncs=sum(temp1);
%K
temp1=[];penk_cpnk=[];
for i=1:6
temp1(i,:)=cell_selecter(Ephys,'drugs',0,'label',1,'geno',7,'sol',1,'optovariant',1,'pair',i);
end
penk_cpnk=sum(temp1);
%% Unlabelled excitatry cell in ntsr1 animal (additional to Penk CCi) only K solution
%CCiK
%K
temp1=[];cci_k=[];
for i=1:6
temp1(i,:)=cell_selecter(Ephys,'drugs',0,'label',0,'geno',6,'sol',1,'optovariant',0,'pair',i);
end
cci_k=sum(temp1);
%K
temp1=[];cci_cpnk=[];
for i=1:6
temp1(i,:)=cell_selecter(Ephys,'drugs',0,'label',1,'geno',6,'sol',1,'optovariant',0,'pair',i);
end
cci_cpnk=sum(temp1);

%% SOM cells 
%CS
temp1=[];som_cs=[];
for i=1:6
temp1(i,:)=cell_selecter(Ephys,'drugs',0,'label',6,'geno',8,'sol',2,'optovariant',1,'pair',i);
end
som_cs=sum(temp1);
%K
temp1=[];som_k=[];
for i=1:6
temp1(i,:)=cell_selecter(Ephys,'drugs',0,'label',6,'geno',8,'sol',1,'optovariant',1,'pair',i);
end
som_k=sum(temp1);

%SOM CPN
%CS
temp1=[];som_cpncs=[];
for i=1:6
temp1(i,:)=cell_selecter(Ephys,'drugs',0,'label',1,'geno',8,'sol',2,'optovariant',1,'pair',i);
end
som_cpncs=sum(temp1);
%K
temp1=[];som_cpnk=[];
for i=1:6
temp1(i,:)=cell_selecter(Ephys,'drugs',0,'label',1,'geno',8,'sol',1,'optovariant',1,'pair',i);
end
som_cpnk=sum(temp1);

%% PV cells 
%CS
temp1=[];pv_cs=[];
for i=1:6
temp1(i,:)=cell_selecter(Ephys,'drugs',0,'label',3,'geno',5,'sol',2,'pair',i);
end
pv_cs=sum(temp1);
%K
temp1=[];pv_k=[];
for i=1:6
temp1(i,:)=cell_selecter(Ephys,'drugs',0,'label',3,'geno',5,'sol',1,'pair',i);
end
pv_k=sum(temp1);

%PV CPN
%CS
temp1=[];pv_cpncs=[];
for i=1:6
temp1(i,:)=cell_selecter(Ephys,'drugs',0,'label',1,'geno',5,'sol',2,'pair',i);
end
pv_cpncs=sum(temp1);
%K
temp1=[];pv_cpnk=[];
for i=1:6
temp1(i,:)=cell_selecter(Ephys,'drugs',0,'label',1,'geno',5,'sol',1,'pair',i);
end
pv_cpnk=sum(temp1);
%remove cell 195 and 200 from ntsr1 cause its pair is with TTX
pv_k([195 200])=0;

%% VIP cells
%CS
temp1=[];vip_cs=[];
for i=1:6
temp1(i,:)=cell_selecter(Ephys,'drugs',0,'label',8,'geno',9,'sol',2,'pair',i);
end
vip_cs=sum(temp1);
%K
temp1=[];vip_k=[];
for i=1:6
temp1(i,:)=cell_selecter(Ephys,'drugs',0,'label',8,'geno',9,'sol',1,'pair',i);
end
vip_k=sum(temp1);

%PV CPN
%CS
temp1=[];vip_cpncs=[];
for i=1:6
temp1(i,:)=cell_selecter(Ephys,'drugs',0,'label',1,'geno',9,'sol',2,'pair',i);
end
vip_cpncs=sum(temp1);
%K
temp1=[];vip_cpnk=[];
for i=1:6
temp1(i,:)=cell_selecter(Ephys,'drugs',0,'label',1,'geno',9,'sol',1,'pair',i);
end
vip_cpnk=sum(temp1);
%%  Read out peaks from train for NTSR CS for train (long), middle frequency (high), highest (highf)
[ntsr_psc_csl] = readout_amp_update(Ephys,ntsr_cs ,1,2,1,2);
[ntsr_psc_cshf] = readout_amp_update(Ephys,ntsr_cs ,2,2,1,2);
[ntsr_psc_cshf2] = readout_amp_update(Ephys,ntsr_cs ,3,2,1,2);
%%  Read out peaks from train for NTSR K for train (long), middle frequency (high), highest (highf)
[ntsr_psc_kl] = readout_amp_update(Ephys,ntsr_k ,1,1,1,2);
[ntsr_psc_khf] = readout_amp_update(Ephys,ntsr_k ,2,1,1,2);
[ntsr_psc_khf2] = readout_amp_update(Ephys,ntsr_k ,3,1,1,2);
%%  Read out peaks from train for CPN NTSR CS for train (long), middle frequency (high), highest (highf)
[ntsr_psc_cpncsl] = readout_amp_update(Ephys,ntsr_cpncs ,1,2,1,2);
[ntsr_psc_cpncshf] = readout_amp_update(Ephys,ntsr_cpncs ,2,2,1,2);
[ntsr_psc_cpncshf2] = readout_amp_update(Ephys,ntsr_cpncs ,3,2,1,2);
%%  Read out peaks from train for CPN NTSR K for train (long), middle frequency (high), highest (highf)
[ntsr_psc_cpnkl] = readout_amp_update(Ephys,ntsr_cpnk ,1,1,1,2);
[ntsr_psc_cpnkhf] = readout_amp_update(Ephys,ntsr_cpnk ,2,1,1,2);
[ntsr_psc_cpnkhf2] = readout_amp_update(Ephys,ntsr_cpnk ,3,1,1,2);

%PENK
%%  Read out peaks from train for PENK CS for train (long), middle frequency (high), highest (highf)
[penk_psc_csl] = readout_amp_update(Ephys,penk_cs ,1,2,1,2);
[penk_psc_cshf] = readout_amp_update(Ephys,penk_cs ,2,2,1,2);
[penk_psc_cshf2] = readout_amp_update(Ephys,penk_cs ,3,2,1,2);
%%  Read out peaks from train for PENK K for train (long), middle frequency (high), highest (highf)
[penk_psc_kl] = readout_amp_update(Ephys,penk_k ,1,1,1,2);
[penk_psc_khf] = readout_amp_update(Ephys,penk_k ,2,1,1,2);
[penk_psc_khf2] = readout_amp_update(Ephys,penk_k ,3,1,1,2);
%%  Read out peaks from train for CPN PENK CS for train (long), middle frequency (high), highest (highf)
[penk_psc_cpncsl] = readout_amp_update(Ephys,penk_cpncs ,1,2,1,2);
[penk_psc_cpncshf] = readout_amp_update(Ephys,penk_cpncs ,2,2,1,2);
[penk_psc_cpncshf2] = readout_amp_update(Ephys,penk_cpncs ,3,2,1,2);
%%  Read out peaks from train for CPN PENK K for train (long), middle frequency (high), highest (highf)
[penk_psc_cpnkl] = readout_amp_update(Ephys,penk_cpnk ,1,1,1,2);
[penk_psc_cpnkhf] = readout_amp_update(Ephys,penk_cpnk ,2,1,1,2);
[penk_psc_cpnkhf2] = readout_amp_update(Ephys,penk_cpnk ,3,1,1,2);
%% CCI
%%  Read out peaks from train for CCI K for train (long), middle frequency (high), highest (highf)
[cci_psc_kl] = readout_amp_update(Ephys,cci_k ,1,1,1,2);
[cci_psc_khf] = readout_amp_update(Ephys,cci_k ,2,1,1,2);
[cci_psc_khf2] = readout_amp_update(Ephys,cci_k ,3,1,1,2);

[cci_psc_cpnkl] = readout_amp_update(Ephys,cci_cpnk ,1,1,1,2);
[cci_psc_cpnkhf] = readout_amp_update(Ephys,cci_cpnk ,2,1,1,2);
[cci_psc_cpnkhf2] = readout_amp_update(Ephys,cci_cpnk ,3,1,1,2);

%SOM
%%  Read out peaks from train for PENK CS for train (long), middle frequency (high), highest (highf)
[som_psc_csl] = readout_amp_update(Ephys,som_cs ,1,2,1,2);
[som_psc_cshf] = readout_amp_update(Ephys,som_cs ,2,2,1,2);
[som_psc_cshf2] = readout_amp_update(Ephys,som_cs ,3,2,1,2);
%%  Read out peaks from train for PENK K for train (long), middle frequency (high), highest (highf)
[som_psc_kl] = readout_amp_update(Ephys,som_k ,1,1,1,2);
[som_psc_khf] = readout_amp_update(Ephys,som_k ,2,1,1,2);
[som_psc_khf2] = readout_amp_update(Ephys,som_k ,3,1,1,2);
%%  Read out peaks from train for CPN PENK CS for train (long), middle frequency (high), highest (highf)
[som_psc_cpncsl] = readout_amp_update(Ephys,som_cpncs ,1,2,1,2);
[som_psc_cpncshf] = readout_amp_update(Ephys,som_cpncs ,2,2,1,2);
[som_psc_cpncshf2] = readout_amp_update(Ephys,som_cpncs ,3,2,1,2);
%%  Read out peaks from train for CPN PENK K for train (long), middle frequency (high), highest (highf)
[som_psc_cpnkl] = readout_amp_update(Ephys,som_cpnk ,1,1,1,2);
[som_psc_cpnkhf] = readout_amp_update(Ephys,som_cpnk ,2,1,1,2);
[som_psc_cpnkhf2] = readout_amp_update(Ephys,som_cpnk ,3,1,1,2);

%PV
%%  Read out peaks from train for PENK CS for train (long), middle frequency (high), highest (highf)
[pv_psc_csl] = readout_amp_update(Ephys,pv_cs ,1,2,1,2);
[pv_psc_cshf] = readout_amp_update(Ephys,pv_cs ,2,2,1,2);
[pv_psc_cshf2] = readout_amp_update(Ephys,pv_cs ,3,2,1,2);
%%  Read out peaks from train for PENK K for train (long), middle frequency (high), highest (highf)
[pv_psc_kl] = readout_amp_update(Ephys,pv_k ,1,1,1,2);
[pv_psc_khf] = readout_amp_update(Ephys,pv_k ,2,1,1,2);
[pv_psc_khf2] = readout_amp_update(Ephys,pv_k ,3,1,1,2);
%%  Read out peaks from train for CPN PENK CS for train (long), middle frequency (high), highest (highf)
[pv_psc_cpncsl] = readout_amp_update(Ephys,pv_cpncs ,1,2,1,2);
[pv_psc_cpncshf] = readout_amp_update(Ephys,pv_cpncs ,2,2,1,2);
[pv_psc_cpncshf2] = readout_amp_update(Ephys,pv_cpncs ,3,2,1,2);
%%  Read out peaks from train for CPN PENK K for train (long), middle frequency (high), highest (highf)
[pv_psc_cpnkl] = readout_amp_update(Ephys,pv_cpnk ,1,1,1,2);
[pv_psc_cpnkhf] = readout_amp_update(Ephys,pv_cpnk ,2,1,1,2);
[pv_psc_cpnkhf2] = readout_amp_update(Ephys,pv_cpnk ,3,1,1,2);
%% VIP
%%  Read out peaks from train for VIP CS for train (long), middle frequency (high), highest (highf)
% [vip_psc_csl] = readout_amp_update(Ephys,vip_cs ,1,2,1,2);
% [vip_psc_cshf] = readout_amp_update(Ephys,vip_cs ,2,2,1,2);
% [vip_psc_cshf2] = readout_amp_update(Ephys,vip_cs ,3,2,1,2);
%%  Read out peaks from train for VIP K for train (long), middle frequency (high), highest (highf)
[vip_psc_kl] = readout_amp_update(Ephys,vip_k ,1,1,1,2);
[vip_psc_khf] = readout_amp_update(Ephys,vip_k ,2,1,1,2);
[vip_psc_khf2] = readout_amp_update(Ephys,vip_k ,3,1,1,2);
%%  Read out peaks from train for CPN VIP CS for train (long), middle frequency (high), highest (highf)
% [vip_psc_cpncsl] = readout_amp_update(Ephys,vip_cpncs ,1,2,1,2);
% [vip_psc_cpncshf] = readout_amp_update(Ephys,vip_cpncs ,2,2,1,2);
% [vip_psc_cpncshf2] = readout_amp_update(Ephys,vip_cpncs ,3,2,1,2);
%%  Read out peaks from train for CPN VIP K for train (long), middle frequency (high), highest (highf)
[vip_psc_cpnkl] = readout_amp_update(Ephys,vip_cpnk ,1,1,1,2);
[vip_psc_cpnkhf] = readout_amp_update(Ephys,vip_cpnk ,2,1,1,2);
[vip_psc_cpnkhf2] = readout_amp_update(Ephys,vip_cpnk ,3,1,1,2);
%% output data
ntsr_cs_cells = [ntsr_psc_csl ntsr_psc_cshf ntsr_psc_cshf2];
ntsr_k_cells = [ntsr_psc_kl ntsr_psc_khf ntsr_psc_khf2];
ntsr_cpncs_cells = [ntsr_psc_cpncsl ntsr_psc_cpncshf ntsr_psc_cpncshf2];
ntsr_cpnk_cells = [ntsr_psc_cpnkl ntsr_psc_cpnkhf ntsr_psc_cpnkhf2];

penk_cs_cells = [penk_psc_csl penk_psc_cshf penk_psc_cshf2];
penk_k_cells = [penk_psc_kl penk_psc_khf penk_psc_khf2];
penk_cpncs_cells = [penk_psc_cpncsl penk_psc_cpncshf penk_psc_cpncshf2];
penk_cpnk_cells = [penk_psc_cpnkl penk_psc_cpnkhf penk_psc_cpnkhf2];


cci_k_cells = [cci_psc_kl cci_psc_khf cci_psc_khf2];
cci_cpnk_cells = [cci_psc_cpnkl cci_psc_cpnkhf cci_psc_cpnkhf2];

som_cs_cells = [som_psc_csl som_psc_cshf som_psc_cshf2];
som_k_cells = [som_psc_kl som_psc_khf som_psc_khf2];
som_cpncs_cells = [som_psc_cpncsl som_psc_cpncshf som_psc_cpncshf2];
som_cpnk_cells = [som_psc_cpnkl som_psc_cpnkhf som_psc_cpnkhf2];

pv_cs_cells = [pv_psc_csl pv_psc_cshf pv_psc_cshf2];
pv_k_cells = [pv_psc_kl pv_psc_khf pv_psc_khf2];
pv_cpncs_cells = [pv_psc_cpncsl pv_psc_cpncshf pv_psc_cpncshf2];
pv_cpnk_cells = [pv_psc_cpnkl pv_psc_cpnkhf pv_psc_cpnkhf2];

%vip_cs_cells = [vip_psc_csl vip_psc_cshf vip_psc_cshf2];
vip_k_cells = [vip_psc_kl vip_psc_khf vip_psc_khf2];
%vip_cpncs_cells = [vip_psc_cpncsl vip_psc_cpncshf vip_psc_cpncshf2];
vip_cpnk_cells = [vip_psc_cpnkl vip_psc_cpnkhf vip_psc_cpnkhf2];

ntsr_k_cells(:,[2 3 5 6 8 9])=ntsr_k_cells(:,[2 3 5 6 8 9])*NaN;
ntsr_cpnk_cells(:,[2 3 5 6 8 9])=ntsr_cpnk_cells(:,[2 3 5 6 8 9])*NaN;

penk_k_cells(:,[2 3 5 6 8 9])=penk_k_cells(:,[2 3 5 6 8 9])*NaN;
penk_cpnk_cells(:,[2 3 5 6 8 9])=penk_cpnk_cells(:,[2 3 5 6 8 9])*NaN;

cci_k_cells(:,[2 3 5 6 8 9])=cci_k_cells(:,[2 3 5 6 8 9])*NaN;
cci_cpnk_cells(:,[2 3 5 6 8 9])=cci_cpnk_cells(:,[2 3 5 6 8 9])*NaN;

som_k_cells(:,[2 3 5 6 8 9])=som_k_cells(:,[2 3 5 6 8 9])*NaN;
som_cpnk_cells(:,[2 3 5 6 8 9])=som_cpnk_cells(:,[2 3 5 6 8 9])*NaN;

pv_k_cells(:,[2 3 5 6 8 9])=pv_k_cells(:,[2 3 5 6 8 9])*NaN;
pv_cpnk_cells(:,[2 3 5 6 8 9])=pv_cpnk_cells(:,[2 3 5 6 8 9])*NaN;

vip_k_cells(:,[2 3 5 6 8 9])=vip_k_cells(:,[2 3 5 6 8 9])*NaN;
vip_cpnk_cells(:,[2 3 5 6 8 9])=vip_cpnk_cells(:,[2 3 5 6 8 9])*NaN;

end
