function [batchopt] = parseExperimentsXls_ephys(path)

[xls_num,xls_txt]=xlsread(path);
%read out colums of interest
animalname    = find(~cellfun(@isempty, strfind(xls_txt(1,:),'AnimalID')));
loadcol        = find(~cellfun(@isempty, strfind(xls_txt(1,:),'BatchAnalyze')));
injectcol        =find(~cellfun(@isempty, strfind(xls_txt(1,:),'InjectionDate')));
datecol       = find(~cellfun(@isempty, strfind(xls_txt(1,:),'ExperimentalDay')));
expcol         = find(~cellfun(@isempty, strfind(xls_txt(1,:),'CellID')));
injecthem= find(~cellfun(@isempty, strfind(xls_txt(1,:),'InjectionH')));
recordhem= find(~cellfun(@isempty, strfind(xls_txt(1,:),'RecordingH')));
genotcol= find(~cellfun(@isempty, strfind(xls_txt(1,:),'Line')));
slicenr= find(~cellfun(@isempty, strfind(xls_txt(1,:),'SliceNr')));
sagflag= find(~cellfun(@isempty, strfind(xls_txt(1,:),'SagFlag')));
label= find(~cellfun(@isempty, strfind(xls_txt(1,:),'Label')));
genotype= find(~cellfun(@isempty, strfind(xls_txt(1,:),'Genotypecode')));
loaddrivecol  = find(~cellfun(@isempty, strfind(xls_txt(1,:),'loaddrive')));

area=find(~cellfun(@isempty, strfind(xls_txt(1,:),'Cortexarea')));
wc=find(~cellfun(@isempty, strfind(xls_txt(1,:),'WC')));
sol=find(~cellfun(@isempty, strfind(xls_txt(1,:),'Solution')));
amp=find(~cellfun(@isempty, strfind(xls_txt(1,:),'Amplifier')));
pair=find(~cellfun(@isempty, strfind(xls_txt(1,:),'Pair')));
disx=find(~cellfun(@isempty, strfind(xls_txt(1,:),'DistX')));
disy=find(~cellfun(@isempty, strfind(xls_txt(1,:),'DistY')));
drug=find(~cellfun(@isempty, strfind(xls_txt(1,:),'Drugs')));
layer=find(~cellfun(@isempty, strfind(xls_txt(1,:),'Layer')));
optovariant=find(~cellfun(@isempty, strfind(xls_txt(1,:),'Optovariant')));
k = 1;

batchopt.XLS.txt = xls_txt;
batchopt.XLS.num = xls_txt;

for i = 2:size(xls_txt,1)
    ana{k}= xls_num(i-1,loadcol-1);
    
    if ~ana{k}
        disp(['skipping experiments ' xls_txt{i,animalname} ' (no batchload flag)']);
        continue
    end
    
    batchopt.mouse{k}     = xls_num(i-1,datecol-1);
    batchopt.mouseID{k}   = xls_txt(i,animalname);
    batchopt.injection{k} = xls_txt(i,injectcol);
    batchopt.line{k} = xls_txt(i,genotcol);
    %batchopt.eye_inj_order{k} = xls_num(i-1,eye_inj_order-1);
    
    expcellids{k}                = xls_txt(i,expcol);
    expcellids2{k}                = xls_txt(i,injecthem);
    expcellids3{k}                = xls_txt(i,recordhem);
    expcellids4{k}                = xls_txt(i,slicenr);
    expcellids5{k}                = xls_txt(i,sagflag);
    expcellids6{k}                = xls_txt(i,label);
    expcellids7{k}                = xls_txt(i,genotype);
    expcellids8{k}                = xls_txt(i,area);
    expcellids9{k}                = xls_txt(i,wc);
    expcellids10{k}                = xls_txt(i,sol);
    expcellids11{k}                = xls_txt(i,amp);
    expcellids12{k}                = xls_txt(i,pair);
    expcellids13{k}                = xls_txt(i,disx);
    expcellids14{k}                = xls_txt(i,disy);
    expcellids15{k}                = xls_txt(i,drug);
    expcellids16{k}                = xls_txt(i,layer);
    expcellids17{k}                = xls_txt(i,optovariant);
    
    batchopt.exp_ids{k}        = str2num((expcellids{k}{1}));
    batchopt.injectH{k}        = (expcellids2{k}{1});
    batchopt.recordH{k}        = (expcellids3{k}{1});
    
    batchopt.slicen{k}         = str2num((expcellids4{k}{1}));
    batchopt.sag{k}            = str2num((expcellids5{k}{1}));
    batchopt.labelcell{k}      = str2num((expcellids6{k}{1}));
    batchopt.geno{k}           = str2num((expcellids7{k}{1}));
    
    batchopt.brainarea{k}           = str2num((expcellids8{k}{1}));
    batchopt.wholecell{k}           = str2num((expcellids9{k}{1}));
    batchopt.solution{k}           = str2num((expcellids10{k}{1}));
    batchopt.amplifier{k}           = str2num((expcellids11{k}{1}));
    batchopt.paired{k}           = str2num((expcellids12{k}{1}));
    batchopt.distancex{k}           = str2num((expcellids13{k}{1}));
    batchopt.distancey{k}           = str2num((expcellids14{k}{1}));
    batchopt.drugs{k}           = str2num((expcellids15{k}{1}));
    batchopt.layer{k}           = str2num((expcellids16{k}{1}));
    batchopt.optovariant{k}           = str2num((expcellids17{k}{1}));
    %batchopt.eye_inj_order{k}  = str2num((expcellids12{k}{1}));
    %     batchopt.spont_ids{k}          = str2num((spontcellids{k}{1}));
    %     batchopt.sftf_ids{k}          = str2num((sftfcellids{k}{1}));
    %     batchopt.puff_ids{k}          = str2num((puffcellids{k}{1}));
    %
    %     batchopt.samesite{k}         = xls_num(i-1,samesitecol-1);
    %     batchopt.baseline{k}         = xls_num(i-1,baselinecol(1)-1);
    %     try
    %         batchopt.baselinepair{k}     = eval(cell2mat(xls_txt(i,basepaircol(1))));
    %         batchopt.baselinepair14{k}   = eval(cell2mat(xls_txt(i,basepair14col)));
    %         batchopt.recovery{k}         = xls_num(i-1,recoverycol-1);
    %     end
    batchopt.loaddrive{k}        = xls_txt(i,loaddrivecol);
    k = k+1;
end

