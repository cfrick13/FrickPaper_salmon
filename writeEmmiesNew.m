
function writeEmmiesNew(director)
% for DD = {'2015_07_10 smFISH'};
CountsString = 'mrna_all.mat';
ConcString = 'mconc.mat';
% for DD = {'2015_01_15 smFISH','2015_01_19 smFISH','2015_01_29 smFISH','2015_01_31 smFISH','2015_03_06 smFISH','2015_03_25 smFISH','2015_03_31 smFISH','2015_04_01 smFISH','2015_08_31 smFISH','2015_09_03 smFISH','2015_12_15 smFISH','2015_12_19 smFISH','2016_01_25 smFISH'};% for DD = {'2015_12_15 smFISH','2015_12_19 smFISH'}
% for DD = {'2016_01_25 smFISH'}
    % for DD = {'2015_01_15 smFISH','2015_01_19 smFISH','2015_01_29 smFISH','2015_01_31 smFISH','2015_03_06 smFISH','2015_03_25 smFISH','2015_03_31 smFISH','2015_04_01 smFISH','2015_05_14 smFISH'};
for DD = {'2015_01_19 smFISH','2015_01_29 smFISH','2015_01_31 smFISH','2015_03_06 smFISH','2015_03_25 smFISH','2015_03_31 smFISH','2015_04_01 smFISH','2015_08_31','2015_09_03','2015_12_15','2015_12_19','2016_01_25','2016_02_09','2016_02_19'}
% for DD = {'2015_01_15 smFISH','2015_01_19 smFISH','2015_01_29 smFISH','2015_01_30 smFISH','2015_01_31 smFISH','2015_03_06 smFISH','2015_03_25 smFISH','2015_03_31 smFISH','2015_04_01 smFISH'}
% for DD = {'2015_01_15','2015_01_19','2015_01_29','2015_01_31','2015_03_06','2015_03_25','2015_03_31','2015_04_01'};    
% for DD = {'2015_08_31','2015_09_03','2015_12_15','2015_12_19','2016_01_25','2016_02_09','2016_02_19'}
% for DD = {'2016_02_09 smFISH','2016_02_19'};

        Dates = char(DD);
Date = Dates(1:10);



A = strcat('D:\Users\zeiss\Pictures\Frick\',Date,' smFISH\');
% CCcellDir = strcat(A,'FISHareaNew7');
% F = 'FISHareaNewMEAN7';
% F = 'FISHareaNewMEAN5';
F = director;
ksize=5;
CCcellDir = strcat(A,F,num2str(ksize(1)));
DICDir = strcat(A,'autoseg\');
COUNTSdir = strcat('D:\Users\zeiss\Documents\MATLAB\AllImagingDataCompiled');
% COUNTSdir = strcat('D:\Users\zeiss\Documents\MATLAB\CountsAndConcManThreshNew');
% CONCdir = strcat('D:\Users\zeiss\Documents\MATLAB\CountsAndConc');
CONCdir = strcat('D:\Users\zeiss\Documents\MATLAB\AllImagingDataCompiled');
mkdir(CONCdir)
cd(CONCdir)

%load the cellLocation info in the CCpositions file
cd (CCcellDir)
filelist = dir(strcat('CCcells*.mat'));
CHANS = findNumberOfVarsInListNormal(filelist,'(594|647)');

cd (DICDir)
filelist = dir(strcat('*sm_*.mat'));
HOURS = findNumberOfVarsInList(filelist,'p[0-9]+');


for CHAN = CHANS
channel = char(CHAN);

for pv = HOURS
cd (DICDir)  
pvalue = char(pv);
disp(pvalue)
filelist = dir(strcat('*sm*',pvalue,'*.mat'));
cfile = char({filelist.name});
CCpositions = load(cfile);

%load 647
cd (CCcellDir)
filelist = dir(strcat('CCcells_',channel,'*',pvalue,'*.mat'));
cfile = char({filelist.name});
EMMIESsavename = cfile;

    if ~isempty(cfile)
    load(cfile); %load CCcells

        for i = 1:length(CCcells.PixelIdxList)
            imgcellz = zeros(CCcells.ImageSize);
            cellz = CCcells.PixelIdxList{i};
            imgcellz(cellz) = 1;
            imgslice = imgcellz(:,:,1);
            CCcells.PixelSlice{i} = find(imgslice>0);
        end

        %channel is already assigned    
        [a,b] = regexp(EMMIESsavename,'(594_snail|594_smad7|594_pai1|647_snail|647_wnt9a|647_ctgf|594_wnt9a|594_ctgf|647_pai1|647_smad7|594_bhlhe40)');
        species = EMMIESsavename(a:b);
        [c,d] = regexp(EMMIESsavename,'(off|low|medlow|med|high|0dot00|0dot02|0dot03|0dot04|0dot07|2dot40)');
        dose = EMMIESsavename(c:d);
        [e,f] = regexp(EMMIESsavename,'[0-9]hr');
        tpoint = EMMIESsavename(e:f);
        %pvalue is already assigned;
        %Date is already assigned;
        % savename= strcat('emmies',chan,'_',Date,'_',tpoint,'_',species,'_',dose,'_',pvalue);

        cellmrnacount = nan(size(CCpositions.cellLocations));
        cellmrnaconcentration =  nan(size(CCpositions.cellLocations));
        cellarea = nan(size(CCpositions.cellLocations));

        cellIDs = CCpositions.cellIDs;
        COUNTS = CCcells.(strcat('mRNA',channel,'inCell'));
        CONC = CCcells.(strcat('mRNA',channel,'inCellconc'));

        for i =1:length(cellIDs)
        cellLocation = CCpositions.cellLocations{i};
        cellID = cellIDs{i};
            if ~isempty(cellLocation)
            ImageSize = CCcells.ImageSize;

                if sum([(cellLocation(1)>2048) (cellLocation(2)>2048)])>0
                    cellLocation =[1 1];
                elseif sum([(cellLocation(1)<0) (cellLocation(2)<0)])>0
                    cellLocation =[1 1];
                else
                    index = sub2ind([ImageSize(1) ImageSize(2)],round(cellLocation(2)),round(cellLocation(1)));
                        for j = 1:CCcells.NumObjects
                            m = CCcells.PixelSlice{j};
                            testtruth = sum(ismember(m,index));
                            if testtruth>0
                                    cellNumber{j} = cellID;
                                    cellmrnacount(i) = COUNTS(j);
                                    cellmrnaconcentration(i) = CONC(j);
                                    cellarea(i) = CCcells.areaOfCell(j);
                            end
                        end

                end
            end
        end



        %save the mRNA counts files
        cd(COUNTSdir)
        if isempty(dir(CountsString)) %write original
            emmiesname = strcat('emmies',channel,'_',Date,'_',tpoint,'_',species,'_',dose,'_',pvalue);
            CCpositions.(emmiesname) = cellmrnacount;
            disp(emmiesname)
            disp(cellmrnacount')
            save(CountsString, '-struct','CCpositions',emmiesname);
        else    %append
            emmiesname = strcat('emmies',channel,'_',Date,'_',tpoint,'_',species,'_',dose,'_',pvalue);
            CCpositions.(emmiesname) = cellmrnacount;
            disp(emmiesname)
            disp(cellmrnacount')
            save(CountsString, '-struct','CCpositions',emmiesname,'-append');
        end

        %save the mRNA concentrations files

        cd(CONCdir)
        if isempty(dir(ConcString))
            emmiesname = strcat('emmies',channel,'_',Date,'_',tpoint,'_',species,'_',dose,'_',pvalue);
            CCpositions.(emmiesname) = cellmrnaconcentration;
            save(ConcString, '-struct','CCpositions',emmiesname);
        else %append
            emmiesname = strcat('emmies',channel,'_',Date,'_',tpoint,'_',species,'_',dose,'_',pvalue);
            CCpositions.(emmiesname) = cellmrnaconcentration;
            save(ConcString, '-struct','CCpositions',emmiesname,'-append');
        end
    end
end
end
end
end






function HOURS = findNumberOfVarsInList(filelist, stringzy)
jjj=1;
for cfile = {filelist.name}
filename = char(cfile);
[aa,bb] = regexp(filename,stringzy);

whostruct = whos('-file',filename);
if length(whostruct)>1
[cc,dd] = regexp(filename,'reference');
refcheck = filename(cc:dd);
if isempty(refcheck)
hours = filename(aa:bb);

if strcmp(hours(1),'s')
    hours(1) = 'p';
end

if jjj==1;
HOURS{jjj} = hours;
jjj=jjj+1;
elseif ~strcmp(HOURS,hours)
    HOURS{jjj} = hours;
    jjj=jjj+1;
else
end
end
end
end
end

function HOURS = findNumberOfVarsInListNormal(filelist, stringzy)
jjj=1;
for cfile = {filelist.name}
filename = char(cfile);
[aa,bb] = regexp(filename,stringzy);

[cc,dd] = regexp(filename,'reference');
refcheck = filename(cc:dd);
if isempty(refcheck)
hours = filename(aa:bb);

if strcmp(hours(1),'s')
    hours(1) = 'p';
end

if jjj==1;
HOURS{jjj} = hours;
jjj=jjj+1;
elseif ~strcmp(HOURS,hours)
    HOURS{jjj} = hours;
    jjj=jjj+1;
else
end
end
end
end
