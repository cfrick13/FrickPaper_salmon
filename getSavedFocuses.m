function getSavedFocuses(Date)
focusStruct=struct();

mfile = mfilename('fullpath');
[~,b] = regexp(mfile,'FrickPaperData');
mfiledir = mfile(1:b+1);
parentdir = mfiledir;
A = parentdir; %directory should end in '...FrickPaperData\' (e.g. D:\FrickPaperData\)

for BB = Date;
% for BB = {'2015_01_19 smFISH','2015_01_29 smFISH','2015_01_31 smFISH','2015_03_06 smFISH','2015_03_25 smFISH','2015_03_31 smFISH','2015_04_01 smFISH'};
ksize = [7];
% kisze = 5;
% for BB = {'2015_01_19 smFISH','2015_01_29 smFISH','2015_01_30 smFISH','2015_01_31 smFISH','2015_03_06 smFISH'};
B = BB{1};
B = strcat(B,' smFISH');

D = '\FLATFIELD';
F = '\FISHareaNew';
G = '\Thresholds';
%make the directories for saving emmies and for saving thresholds
mkdir (strcat(A,B,F,num2str(ksize(1))));
mkdir (strcat(A,B,G,num2str(ksize(1))));
THR = (strcat(A,B,G,num2str(ksize(1))));
EXP = strcat(A,B,D); 
cd (EXP)

%determine number of channels used
filelist = dir(strcat('*','*.tif'));
CHANSS = findNumberOfVarsInList(filelist,'(Fluor 594|Fluor 647)');
CHANS = cellfun(@(x) x(end-2:end),CHANSS,'UniformOutput',0);

%determine number of pvalues
expfilelist = dir(strcat('*DIC*z10*.tif'));
expfilenames = {expfilelist.name};
PVALUES = findNumberOfVarsInListP(expfilelist,'p[0-9]+');


% PVALUES = {'p24'};
i=1;
for pvaluecell = PVALUES
pvalue = char(pvaluecell);


[Details,Fishdotstack,Imgstack] = makeimagestack(EXP,pvalue,ksize);

%determine what channels are present
Fnames = fieldnames(Details);
fname=Fnames{1};
[aa,bb] = regexp(fname,'(594|647)');
if ~isempty(aa)
fishdotstack = Fishdotstack.(char(fname));
imgstack = Imgstack.(char(fname));
details = Details.(char(fname));    
    
%Define the minimum size of dots that qualify as mRNA
mrnasize = 6;

cd(THR)
focusthresh = 12000;
onlysegmented = fishdotstack>(focusthresh);
dotss = zeros(1,size(fishdotstack,3));
parfor ii=1:size(fishdotstack,3)
% for ii=1:size(fishdotstack,3)
onlysegmentedslice = onlysegmented(:,:,ii);
mmm = bwconncomp(onlysegmentedslice,8);
pix = mmm.PixelIdxList;
pixies = cellfun(@length,pix);
pixx = pixies>mrnasize;
dotss(ii) = sum(pixx);
end



focusPoint = find(dotss==max(dotss),1,'first'); %determine the most in focus frame
focusStruct(i).stack = dotss;
disp(focusStruct)
focusStruct(i).pvalue = pvalue;
focusStruct(i).peak = focusPoint;
i=i+1;


save(strcat('focus_',pvalue,'_',details.species,'.mat'),'focusPoint');


end
end
cd('D:\Users\zeiss\Documents\MATLAB\')
save(strcat(B(1:10),'focusStruct.mat'),'focusStruct');
focusStruct=struct();
stophere=1;
end
end

function [Details,Fishdotstack,Stack] = makeimagestack(EXP,pvalue,ksize)
cd (EXP)
Stack = struct();
Fishdotstack = struct();

for chan = {'Fluor 594'}
imagefilelist = dir(strcat('*',pvalue,'-*',char(chan),'*.tif'));
if ~isempty(imagefilelist)
imagefilenames = {imagefilelist.name};
imageinfo = imfinfo(char(imagefilenames{1}));
imgx = imageinfo.Width;
imgy = imageinfo.Height;
imgz = length({imagefilelist.name});
dims = [imgx imgy imgz];
details.dims = dims;
stack = zeros(dims(1),dims(2),dims(3));
i=1;
for cfile = {imagefilelist.name}
    stack(:,:,i) = imread(char(cfile));
    i=i+1;
end
fishdotstack = LaplacianOfGaussianStack(stack,dims,ksize); %function in code




filename = char(imagefilenames{1});
[aa,bb] = regexp(filename,'p[0-9]+');
details.pvalue = filename(aa:bb);
[cc,dd] = regexp(filename,'[0-9]hr');
details.tpoint = filename(cc:dd);
[ee,ff] = regexp(filename,'(off|low|medlow|med|high|0dot00|0dot02|0dot03|0dot04|0dot07|2dot40)');
dose = filename(ee:ff);
if strcmp(chan,'Fluor 594')
[gg,hh] = regexp(filename,'(594_snail|594_smad7|594_pai1|594_pmepa1|594_tieg|594_bhlhe40|594_wnt9a|594_ctgf)');
species = filename(gg:hh);
    if isempty(species)
    [mm,nn] = regexp(filename,'(smad7|snail|row1|row2|row3)');
    species = filename(mm:nn);
    end
elseif strcmp(chan,'Fluor 647')
[gg,hh] = regexp(filename,'(647_snail|647_smad7|647_pai1|647_wnt9a|647_ctgf)');    
species = filename(gg:hh);
    if isempty(species)
    species = 'pai1';
    end
end



details.chanspecie = species;
if strcmp(details.chanspecie,'row1')||strcmp(details.chanspecie,'row2')||strcmp(details.chanspecie,'row3')
    details.chanspecie = 'snail';
end

[yy,zz] = regexp(filename,'(row1|row2|row3)');
row = filename(yy:zz);
    
[ii,jj] = regexp(filename,'flat_2');
details.date = filename(jj:jj+9);

details.channelnumber = chan;

if strcmp(details.tpoint,'0hr')
    dose = '0dot00';
end
% 
% %%%dose naming correction
% if strcmp(dose,'off')
%     dose = '0dot00';
% elseif strcmp(dose,'low')
%     dose = '0dot02';
% elseif strcmp(dose,'medlow')
%     dose = '0dot04';
% elseif strcmp(dose,'med')
%     dose = '0dot07';
% elseif strcmp(dose,'high')
%     dose = '2dot40';
% end

if strcmp(details.date,'2015_1_30 ')
    details.date = '2015_01_30';
%%%dose correction for 2015_01_30
if strcmp(chan,'Fluor 647')
    
    if strcmp(row,'row1')
    dose = '2dot40';    
    elseif strcmp(row,'row2')
    dose = '0dot07';
    end
    species = '647_pai1';
    
elseif strcmp(chan,'Fluor 594')
    
    if strcmp(species,'row1')
    dose = '2dot40';
    elseif strcmp(species,'row2')
    dose = '0dot07';
    end
    species = '594_snail';
end

if strcmp(details.tpoint,'0hr')
    dose = '0dot00';
end

end


if strcmp(details.date,'2014_12_2 ')
    if strcmp(species,'row2')
        species = '594_snail';
    elseif strcmp(species,'row3')
        species = '594_snail';
    end
end
        
    
       
%%%species naming correction
if strcmp(chan,'Fluor 647')
    if strcmp(species,'smad7')
    species = '647_pai1';
    elseif strcmp(species,'snail')
    species = '647_pai1';
    elseif strcmp(species,'pai1')
    species = '647_pai1';
    end
else
    if strcmp(species,'smad7')
    species = '594_smad7';
    elseif strcmp(species,'snail')
    species = '594_snail';
    elseif strcmp(species,'pai1')
    species = '594_pai1';
    
    elseif strcmp(species,'pmepa1')
        species = '594_pmepai1';
    elseif strcmp(species,'tieg')
        species = '594_tieg';
    elseif strcmp(species,'bhlhe40')
        species = '594_bhlhe40';
    end
end

details.dose = dose;
details.species = species;

disp(details)

dims = details.dims;


Stack.(strcat('i',details.species,'i'))=stack;
Fishdotstack.(strcat('i',details.species,'i'))=fishdotstack;
Details.(strcat('i',details.species,'i')) = details;
end
end
end

function  LoGstack = LaplacianOfGaussianStack(imgstack,dims,ksize)
        LoGstack = zeros(dims(1),dims(2),dims(3));
        for i = 1:size(imgstack,3)
        LoGstack(:,:,i) = logMasked(imgstack(:,:,i),ksize);
        end
        
end

function kernel = chooseKernel(ksize)
if ksize ==5
kernel = [-4 -1  0 -1 -4;...
     -1  2  3  2 -1;...
      0  3  4  3  0;...
     -1  2  3  2 -1;...
     -4 -1  0 -1 -4];


% % % -4 -1  0 -1 -4
% % % -1  2  3  2 -1
% % % 0  3  4  3  0
% % % -1  2  3  2 -1
% % % -4 -1  0 -1 -4

elseif ksize == 7
kernel =[-10 -5 -2 -1 -2 -5 -10;... 
    -5  0  3  4  3  0  -5;... 
    -2  3  6  7  6  3  -2;... 
    -1  4  7  8  7  4  -1;... 
    -2  3  6  7  6  3  -2;... 
    -5  0  3  4  3  0  -5;... 
    -10 -5 -2 -1 -2 -5 -10];... 
    
% % % % -10 -5 -2 -1 -2 -5 -10 
% % % % -5  0  3  4  3  0  -5 
% % % % -2  3  6  7  6  3  -2 
% % % % -1  4  7  8  7  4  -1 
% % % % -2  3  6  7  6  3  -2 
% % % % -5  0  3  4  3  0  -5
% % % % -10 -5 -2 -1 -2 -5 -10
end
end

function bw = logMasked(im,ksize,varargin)
% Discrete Laplacian
kernel = chooseKernel(ksize);
%% image filtering
lapFrame = imfilter(im,kernel,'repl');
if ~isempty(varargin)
    bw=lapFrame.*uint16(varargin{1}>0);
else
    bw=lapFrame;
end
end

function HOURS = findNumberOfVarsInList(filelist, stringzy)
jjj=1;
HOURS=[];
for cfile = {filelist.name}
filename = char(cfile);
[aa,bb] = regexp(filename,stringzy);
hours = filename(aa:bb);
[cc,dd] = regexp(filename,'reference');
ref = filename(cc:dd);

if ~isempty(hours)
if isempty(ref)
if jjj==1;
HOURS{jjj} = hours;
jjj=jjj+1;
elseif ~strcmp(HOURS,hours)
    HOURS{jjj} = hours;
    jjj=jjj+1;
else
end

else
end
end
end
end


function HOURS = findNumberOfVarsInListP(filelist, stringzy)
jjj=1;
HOURS=[];
for cfile = {filelist.name}
filename = char(cfile);
[aa,bb] = regexp(filename,stringzy);
hours = filename(aa:bb);
[cc,dd] = regexp(filename,'reference');
ref = filename(cc:dd);
if isempty(ref)
HOURS{jjj} = hours;
jjj=jjj+1;
end
end
HOURS = unique(HOURS);
end

