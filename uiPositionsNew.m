function uiPositionsNew
global fff CellImage correct DICimage AnnotationsJPGDir pvalueold stats A Position mm CCcellDir DICDir CCcells ExcelDir cellNumbers HOURS xspot AnnotationsDir ff f
close all
pvalueold ='old';
correct = 0;
Date = '2016_02_09';
% Date = '2015_01_31';

stats=struct();
DICimage = zeros(2048,2048);
basefolder = 'F:\FrickPaper\';
A = strcat(basefolder,Date,' smFISH\');
ExcelDir = strcat(basefolder,'\ImagingResults');
CCcellDir = strcat(A,'FISHareaNew7');
DICDir = strcat(A,'autoseg\');
AnnotationsDir = strcat(A,'ANNOTATIONS\');
AnnotationsJPGDir = strcat(A,'ANNOTATIONSJPG\');

%allow for opening all segmented images
cd(DICDir)
primarylist = dir('*.tif');
HOURS = findNumberOfVarsInList(primarylist, 'p[0-9]+');

% %only allow to open images that have excel data from movies
% cd(ExcelDir)
% date = A(end-17:end-8);
% filelist = dir(strcat(date,'*SHEET.xlsx'));
% filename = char({filelist.name});
% 
% [status,sheets] = xlsfinfo(filename)
% HOURS =  findNumberOfVarsInListExcel(sheets, 's[0-9]+');


ff= figure;
ff.Visible ='off';
ff.Position =[1100,10,1000,1000];

fff = figure;
fff.Visible ='off';
fff.Position = [1000,100,1000,1000];

f = figure;
f.Visible ='off';
f.Position =[10,10,800,800];
xspot = 670;

htext = uicontrol('Style','text','String','Choose Position to Correct',...
          'Position',[xspot-50,750,200,15]);
hpopup = uicontrol('Style','popupmenu',...
          'String',HOURS',...
          'Position',[xspot,720,100,25],...
          'Callback',@popup_menu_Callback);
htext = uicontrol('Style','text','String','Choose Cell# to Identify',...
          'Position',[xspot-50,680,200,15]);      

% hNext = uicontrol('Style','pushbutton','String','Next',...
%           'Position',[xspot,600,70,25],...
%           'Callback',@Nextbutton_Callback);
hSetPos = uicontrol('Style','pushbutton','String','SetPositionFromAnnotation',...
          'Position',[xspot+50,650,70,25],...
          'Callback',@SetPosbutton_Callback);
hSavePositionInfo = uicontrol('Style','pushbutton','String','SavePositionInfo',...
          'Position',[xspot,300,70,25],...
          'Callback',@SavePositionInfo_Callback);
hupdatePositions = uicontrol('Style','pushbutton','String','UpdatePositions',...
          'Position',[xspot-40,500,70+40,75],...
          'Callback',@updatePostions_Callback);      

%    ha = axes(ax,'Units','Pixels','Position',[50,60,200,185]); 
align([hpopup,htext,],'Center','None');
      
      

      
f.Visible = 'on'   ;
CellImage = zeros(2048,2048);
   f.Units = 'normalized';
   htext.Units = 'normalized';
   hpopup.Units = 'normalized';
   hSetPos.Units = 'normalized';
   hSplit.Units = 'normalized';
   hSaveImage.Units = 'normalized';
   hAdd.Units = 'normalized';
   hRemoveArea.Units = 'normalized';

% ha = axes('Units','Pixels','Position',[50,60,200,185]); 
% mm = subplot(1,2,1);
mm = axes;

Position = [0.05 0.05 0.7 0.7];
mm.Position = Position;
himg = imshow(CellImage);

set(f,'KeyPressFcn',@keypress);
end

function keypress(fig_obj,~)
global  ImageDetails
key = get(fig_obj,'CurrentKey');
    switch key
        case 'q'
            Previousbutton_Callback([],[])
        case 'w'
            Nextbutton_Callback([],[])
        case 's'
            SavePositionInfo_Callback([],[])
        case 'v'
            updatePostions_Callback([],[])
    end
end
function updatePostions_Callback(~,~)
global correct
correct=1;
loadCellInfoAndImages
correct=0;
end


function popup_menu_Callback(source,~)
global pvalue cellNumbers xspot f correct pvalueold DICDir
% Determine the selected data set.
 str = source.String;
 val = source.Value;
 pvalue = char(str{val});


correct=0;    

    
    
loadCellInfoAndImages
hpopupCellID = uicontrol('Style','popupmenu',...
          'String',cellNumbers',...
          'Position',[xspot-20,100,40,25],...
          'Callback',@popup_menuCellID_Callback);
      
hpopupCellIDTwo = uicontrol('Style','popupmenu',...
          'String',cellNumbers',...
          'Position',[xspot-20,650,45,30],...
          'Callback',@popup_menuCellIDTwo_Callback);
end


function popup_menuCellID_Callback(source,~)
global cellIDnumber CCcells CellImage f val
% Determine the selected data set.
 str = source.String;
 val = source.Value;
 cellIDnumber = char(str{val});
 
 SetPosbutton_Callback(1,1) 


end
function SetPosbutton_Callback(~,~)
 global  CellImage cellIDnumber f CCpositions val correct
 %numberOfCells
 %directory with movie spreadsheets
 %click to label cell corresponding to number
 %label the cell with that number
 
 correct=0;
  cellIDs = CCpositions.cellIDs;
  cellLocations = CCpositions.cellLocations;
 
   % Display mesh plot of the currently selected data.
     [cellx,celly] = ginput(1);
      output.cellx=cellx;
      output.celly=celly;
       

    cellIDs{val} = num2str(cellIDnumber);
    cellLocations{val} = [cellx celly];

    CCpositions.cellIDs = cellIDs;
    CCpositions.cellLocations = cellLocations;
    
    updateImage
end


function popup_menuCellIDTwo_Callback(source,~)
global cellIDnumber CCcells CellImage f val
% Determine the selected data set.
 str = source.String;
 val = source.Value;
 cellIDnumber = char(str{val});
 
 SetPosbuttonTwo_Callback(1,1) 
end
function SetPosbuttonTwo_Callback(~,~)
 global  CellImage cellIDnumber f CCpositions val correct
 %numberOfCells
 %directory with movie spreadsheets
 %click to label cell corresponding to number
 %label the cell with that number
 
 correct=0;
  cellIDs = CCpositions.cellIDs;
  cellLocations = CCpositions.cellLocations;
 
   % Display mesh plot of the currently selected data.
     [cellx,celly] = ginput(1);
      output.cellx=cellx;
      output.celly=celly;
       

    cellIDs{val} = num2str(cellIDnumber);
    cellLocations{val} = [cellx celly];

    CCpositions.cellIDs = cellIDs;
    CCpositions.cellLocations = cellLocations;
    
    updateImage
    SavePositionInfo_Callback([],[])
    correct =1;
    loadCellInfoAndImagesTwo
    SavePositionInfo_Callback([],[])
end


function loadCellInfoAndImagesTwo
global CellImage correct DICimage cellIDnumber savename CCcellDir DICDir pvalue ExcelDir cellNumbers A HOURS f CCpositions

%load excel file to get info about which cells are present
cd(ExcelDir)
date = A(end-17:end-8);filelist = dir(strcat(date,'*SHEET.xlsx'));filename = char({filelist.name});
svalue = pvalue;svalue(1) = 's';

    if ismember(pvalue,HOURS)
        [spreadsheet,~,~] = xlsread(filename,svalue);spreadBacks = spreadsheet(3:14,:);
        kkk = find(spreadBacks<10 & spreadBacks>0);
        [~,zz] = ind2sub(size(spreadBacks),kkk);
            if isempty(zz)
                endhere=1;
            end
        spreadTraces = spreadsheet(10:end,:);jjj = find(spreadTraces<10 & spreadTraces>0);[a,b] = ind2sub(size(spreadTraces),jjj);cellNumberz = spreadsheet(a(1)+9:2:end,b(1));
        cellNumbers = cell(length(cellNumberz),1);
            for i = 1:length(cellNumberz)
                cellNumbers{i} = num2str(cellNumberz(i));
            end
    else
        prompt = cell(1,1);
        prompt{1} = 'numbers of cells?';
        Ins = inputdlg(prompt,strcat('enter numbers separated by a space'),[1 50]);

        cellNumString = Ins{1};
        cellNumberz = str2double(cellNumString);
        cellNumbers = strsplit(cellNumString);
    end



cd(DICDir)
%load autoseg corrected segmentation .mats
filelist = dir(strcat('*sm_',pvalue,'*.mat'));
cfile = char({filelist.name});
CCpositions = load(cfile); %this loads segLogical and potentially the cell positions information, if present
fnames = fieldnames(CCpositions);
    if ismember('cellLocations',fnames)
        if correct == 1
            xynew = vertcat(CCpositions.cellLocations{:});    
            scalefactor = 4.9;
            cd ..
            cd('CENTROIDS')
            svalue = pvalue; svalue(1)='s';
            load(strcat('CENTROIDS-',svalue,'.mat'));
            xy = CENTROIDS.centroids;
            xyscaled=xy.*scalefactor;

            %%%%%%%%%%%%%%%%%%%
%             whichcell = str2num(char(inputdlg('which cell are you locating?','choose cell to correct locations',1)));
            whichcell = str2num(cellIDnumber)
%             [x,y] = ginput(1);
            cellLocations = CCpositions.cellLocations;
            xy = cellLocations{whichcell};
            x=xy(1);
            y=xy(2);
            %%%%%%%%%%%%%%%%%%%%
            xadj = xyscaled(whichcell,1)-x;
            yadj = xyscaled(whichcell,2)-y;
            xyadj(:,1) = xyscaled(:,1)-xadj;
            xyadj(:,2) = xyscaled(:,2)-yadj;

            [Idx,Eps] = knnsearch(xynew,xyadj,'K',3);
            EpsOne = Eps(:,1);
            belowthresh = (EpsOne<200);
            belowthreshidx = (belowthresh==0);
%             xynewch = xynew;
            xynewch = xyadj;
            xynewch(belowthresh,1) = xynew(Idx(belowthresh),1);
            xynewch(belowthresh,2) = xynew(Idx(belowthresh),2);
            xynewch(belowthreshidx,1) = xyadj(belowthreshidx,1);
            xynewch(belowthreshidx,2) = xyadj(belowthreshidx,2);
            xxyy = cell(1,size(xynewch,1));
                for i=1:length(xxyy)
                   xxyy(i) = {xynewch(i,:)};
                   cellIDs{i} = num2str(i);
                end
            CCpositions.cellLocations = xxyy;
            CCpositions.cellIDs = cellIDs;
        end
    else
        cellLocations = cell(length(cellNumberz),1);
        cellIDs = cell(length(cellNumberz),1);
        CCpositions.cellLocations = cellLocations;  
        CCpositions.cellIDs = cellIDs;
    end


%make the label matrix of image
CCimage = bwconncomp(CCpositions.segLogical);
CellImage = labelmatrix(CCimage);
savename = cfile;

%load DIC image
cd(DICDir)
    DICfilelist = dir(strcat('*sm_DIC_',pvalue,'*'));
    DICimage = imread(char(DICfilelist.name));
cd(CCcellDir)
    updateImage
end
function loadCellInfoAndImages
global CellImage correct DICimage savename CCcellDir DICDir pvalue ExcelDir cellNumbers A HOURS f CCpositions

%load excel file to get info about which cells are present
cd(ExcelDir)
date = A(end-17:end-8);filelist = dir(strcat(date,'*SHEET.xlsx'));filename = char({filelist.name});
svalue = pvalue;svalue(1) = 's';

    if ismember(pvalue,HOURS)
        [spreadsheet,~,~] = xlsread(filename,svalue);spreadBacks = spreadsheet(3:14,:);
        kkk = find(spreadBacks<10 & spreadBacks>0);
        [~,zz] = ind2sub(size(spreadBacks),kkk);
            if isempty(zz)
                endhere=1;
            end
        spreadTraces = spreadsheet(10:end,:);jjj = find(spreadTraces<10 & spreadTraces>0);[a,b] = ind2sub(size(spreadTraces),jjj);cellNumberz = spreadsheet(a(1)+9:2:end,b(1));
        cellNumbers = cell(length(cellNumberz),1);
            for i = 1:length(cellNumberz)
                cellNumbers{i} = num2str(cellNumberz(i));
            end
    else
        prompt = cell(1,1);
        prompt{1} = 'numbers of cells?';
        Ins = inputdlg(prompt,strcat('enter numbers separated by a space'),[1 50]);

        cellNumString = Ins{1};
        cellNumberz = str2double(cellNumString);
        cellNumbers = strsplit(cellNumString);
    end



cd(DICDir)
%load autoseg corrected segmentation .mats
filelist = dir(strcat('*sm_',pvalue,'*.mat'));
cfile = char({filelist.name});
CCpositions = load(cfile); %this loads segLogical and potentially the cell positions information, if present
fnames = fieldnames(CCpositions);
    if ismember('cellLocations',fnames)
        if correct == 1
            xynew = vertcat(CCpositions.cellLocations{:});    
            scalefactor = 4.9;
            cd ..
            cd('CENTROIDS')
            svalue = pvalue; svalue(1)='s';
            load(strcat('CENTROIDS-',svalue,'.mat'));
            xy = CENTROIDS.centroids;
            xyscaled=xy.*scalefactor;

            %%%%%%%%%%%%%%%%%%%
            whichcell = str2num(char(inputdlg('which cell are you locating?','choose cell to correct locations',1)));
            [x,y] = ginput(1);
            %%%%%%%%%%%%%%%%%%%%
            xadj = xyscaled(whichcell,1)-x;
            yadj = xyscaled(whichcell,2)-y;
            xyadj(:,1) = xyscaled(:,1)-xadj;
            xyadj(:,2) = xyscaled(:,2)-yadj;

            [Idx,Eps] = knnsearch(xynew,xyadj,'K',3);
            EpsOne = Eps(:,1);
            belowthresh = (EpsOne<200);
            belowthreshidx = (belowthresh==0);
%             xynewch = xynew;
            xynewch = xyadj;
            xynewch(belowthresh,1) = xynew(Idx(belowthresh),1);
            xynewch(belowthresh,2) = xynew(Idx(belowthresh),2);
            xynewch(belowthreshidx,1) = xyadj(belowthreshidx,1);
            xynewch(belowthreshidx,2) = xyadj(belowthreshidx,2);
            xxyy = cell(1,size(xynewch,1));
                for i=1:size(xxyy,2)
                   xxyy(i) = {xynewch(i,:)};
                   cellIDs{i} = num2str(i);
                end
            CCpositions.cellLocations = xxyy;
            CCpositions.cellIDs = cellIDs;
        end
    else
        cellLocations = cell(length(cellNumberz),1);
        cellIDs = cell(length(cellNumberz),1);
        CCpositions.cellLocations = cellLocations;  
        CCpositions.cellIDs = cellIDs;
    end


%make the label matrix of image
CCimage = bwconncomp(CCpositions.segLogical);
CellImage = labelmatrix(CCimage);
savename = cfile;

%load DIC image
cd(DICDir)
    DICfilelist = dir(strcat('*sm_DIC_',pvalue,'*'));
    DICimage = imread(char(DICfilelist.name));
cd(CCcellDir)
    updateImage
end
function updateImage
global fff CellImage DICimage savename Position mm AnnotationsDir pvalue ff f CCpositions AnnotationsJPGDir
% CC = bwconncomp(CellImage);
% CellImage = labelmatrix(CC);
% stats = regionprops(CellImage,'PixelIdxList');

% for i=1:length(stats)
%     area= length(stats(i).PixelIdxList);
%     m = stats(i).PixelIdxList;
%     if area>10000
%     else
%         CellImage(m) = 0;
%     end
% end
figure(ff)
axes(mm);
children = findobj(mm,'Type','image');
delete(children);
children = findobj(mm,'Type','text');
delete(children);
mm.Position = Position;
mm.NextPlot = 'replace';
rgb = label2rgb(CellImage,'jet','k','shuffle');hold on
imshow(rgb);
title(savename);
j = imagesc(DICimage);
j.AlphaData = 0.8;
colormap('gray')
h = gca;
h.CLim=([median(median(DICimage))./2 median(median(DICimage)).*2]);


for i=1:length(CCpositions.cellLocations)
    centroid = CCpositions.cellLocations{i};
    cellID = CCpositions.cellIDs{i};
    
    if ~isempty(centroid)
    text(centroid(1)+1,centroid(2)+1,char(cellID),'FontSize',20,'Color',[0 0 0]);hold on
    text(centroid(1),centroid(2),char(cellID),'FontSize',18,'Color',[1 1 1]);hold on
    end
end
% end

mm.NextPlot = 'replace';


figure(ff)
cd (AnnotationsDir)
cd ..
    if isempty(dir('*ANNOTATIONSJPG'));
        scene = 'scene*00';
        scene(end-1:end) = pvalue(end-1:end);
        filelist = dir(strcat('*',scene,'*.tif'));
        
        annfilenames = {filelist.name};
        if ~isempty(annfilenames)
        ann = loadUpFinalImageOfStack(char(annfilenames{1}));
        else
            cd('ANNOTATIONS')
        scene = 's00';
        scene(end-1:end) = pvalue(end-1:end);
        filelist = dir(strcat('*',scene,'*.jpg'));
        annfilenames = {filelist.name};
        ann = imread(char(annfilenames));
        cd .. 
        end
            
    else
        cd (AnnotationsJPGDir)
        svalue = pvalue;svalue(1) = 's';
        filelist = dir(strcat('*',svalue,'*.jpg'));
        annfilenames = {filelist.name};
        ann = imread(char(annfilenames));
    end

ff.Visible = 'on';
imshow(ann)
imgff = gca;
imgff.YDir = 'reverse';

% figure(fff)
cd (AnnotationsDir)
cd ..
%     if ~isempty(dir('*DOTTED*'));
%         cd('DOTTED')
% %         fff.Visible = 'on';
%         ff.Visible = 'on';
%         scene = 'p00';
%         scene(end-1:end) = pvalue(end-1:end);
%         filelist = dir(strcat('*DOTTED*',scene,'*'));
%         annfilenames = {filelist.name};
%         DotStack = loadUpImageStack(annfilenames);
%         zproj = max(DotStack,[],3);
%         imgfff=gca;
%         imgfff.YDir = 'reverse';
%         
%         figure(ff)
%         imshow(rgb);
%         hold on;
%         j = imagesc(logical(zproj));
%         j.AlphaData = 0.8;
%         colormap('gray')
%         h = gca;
% 
%     end
%     
    cd (AnnotationsDir)
cd ..
    
       


figure(f)

end



function Stack=loadUpImageStack(filenames)
cfile = filenames;
FileTif = char(cfile);
InfoImage=imfinfo(FileTif);

mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
NumberImages=length(InfoImage);
Stack=zeros(nImage,mImage,NumberImages,'uint16');
 
TifLink = Tiff(FileTif, 'r');
bb=1;
i = 1;

for i = 1:NumberImages;

   TifLink.setDirectory(i);
   FinalImage=TifLink.read();
   Stack(:,:,i) = FinalImage;
end
 
end
function FinalImage=loadUpFinalImageOfStack(filenames)
cfile = filenames;
FileTif = char(cfile);
InfoImage=imfinfo(FileTif);

mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
NumberImages=length(InfoImage);
FinalImage=zeros(nImage,mImage,NumberImages,'uint16');
 
TifLink = Tiff(FileTif, 'r');
bb=1;
i = 1;

i = NumberImages;

   TifLink.setDirectory(i);
   FinalImage=TifLink.read();
   
 
end


function  Nextbutton_Callback(~,~)
global  pvalue correct cellNumbers xspot

pnum = str2double(pvalue(2:3));
pnum = pnum+1;
pstr = num2str(pnum);
if length(pstr)==1
    pvalue = 'p00';
    pvalue(3) = pstr;
else
    pvalue = 'p00';
    pvalue(2:3) = pstr;
end

correct=0;
loadCellInfoAndImages
hpopupCellID = uicontrol('Style','popupmenu',...
          'String',cellNumbers',...
          'Position',[xspot-20,100,40,25],...
          'Callback',@popup_menuCellID_Callback);
      
hpopupCellIDTwo = uicontrol('Style','popupmenu',...
          'String',cellNumbers',...
          'Position',[xspot-20,650,45,30],...
          'Callback',@popup_menuCellIDTwo_Callback);
end
function  Previousbutton_Callback(~,~)
global  pvalue correct cellNumbers xspot

pnum = str2double(pvalue(2:3));
pnum = pnum-1;
pstr = num2str(pnum);
if length(pstr)==1
    pvalue = 'p00';
    pvalue(3) = pstr;
else
    pvalue = 'p00';
    pvalue(2:3) = pstr;
end

correct=0;
loadCellInfoAndImages
hpopupCellID = uicontrol('Style','popupmenu',...
          'String',cellNumbers',...
          'Position',[xspot-20,100,40,25],...
          'Callback',@popup_menuCellID_Callback);
      
hpopupCellIDTwo = uicontrol('Style','popupmenu',...
          'String',cellNumbers',...
          'Position',[xspot-20,650,45,30],...
          'Callback',@popup_menuCellIDTwo_Callback);
end


function SavePositionInfo_Callback(~,~)
global A pvalue CCcells DICDir savename CCpositions
cd(DICDir)

cellLocations = CCpositions.cellLocations;
cellIDs = CCpositions.cellIDs;

save(strcat(savename),'cellLocations','cellIDs','-append');
disp(strcat('saved...','cellLocations and cellIDs_',pvalue,'.mat'))
% savethatimage(savename,CellImage,A)
end
function savethatimage(savename,CellImage,A)
cd (strcat(A,'autoseg\'));
imwrite(CellImage,strcat(savename),'Tiff');

segLogical = logical(CellImage);
save(strcat(savename(1:end-4),'.mat'),'segLogical');

disp(strcat('image ... ',savename,'...saved'));
end

function HOURS = findNumberOfVarsInList(filelist, stringzy)
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