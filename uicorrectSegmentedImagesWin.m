function uicorrectSegmentedImagesWin
global If DICimage stats A Position mm AnnotationsDir ff f SegLocation HOURS ffposition
% close all
%%%%%%%%%%%%%%%%%%%%%%%%
Date = '2016_02_09';
% Date = '2016_02_09';
%%%%%%%%%%%%%%%%%%%%%%%%

stats=struct();
DICimage = zeros(2048,2048);

mfile = mfilename('fullpath');
[a,b] = regexp(mfile,'FrickPaperData');
mfiledir = mfile(1:b+1);
parentdir = mfiledir;
A = strcat(parentdir,Date,' smFISH\');
    % A = strcat('D:\Users\zeiss\Pictures\Frick\',Date,' smFISH\');
    % A = strcat('F:\FrickPaper\',Date,' smFISH\');
    % A = strcat('\Users\frick\Documents\Goentoro_Lab\DATA\current\',Date,' smFISH\');
% AnnotationsDir = strcat(A,'ANNOTATIONS\');
AnnotationsDir = strcat(A,'ANNOTATIONS');
SegLocation = (strcat(A,'autoseg\'));
cd(SegLocation)

primarylist = dir('*.tif');
HOURS = findNumberOfVarsInList(primarylist, 'p[0-9]+');

f = figure;
f.Visible ='off';
f.Position =[10,10,1200,1100];
xspot = 1050;

htext = uicontrol('Style','text','String','Choose Position to Correct',...
          'Position',[xspot,550,60,15]);
hpopup = uicontrol('Style','popupmenu',...
          'String',HOURS',...
          'Position',[xspot,500,100,25],...
          'Callback',@popup_menu_Callback);
hAdd = uicontrol('Style','pushbutton','String','AddArea',...
          'Position',[xspot,600,70,25],...
          'Callback',@addareabutton_Callback);
hRemoveArea = uicontrol('Style','pushbutton','String','RemoveArea',...
          'Position',[xspot,650,70,25],...
          'Callback',@removeareabutton_Callback);
hDelete = uicontrol('Style','pushbutton','String','Delete',...
          'Position',[xspot,700,70,25],...
          'Callback',@deletebutton_Callback);
hSplit = uicontrol('Style','pushbutton',...
          'String','Split',...
          'Position',[xspot,750,70,25],...
          'Callback',@splitbutton_Callback);
hFillHoles = uicontrol('Style','pushbutton',...
          'String','FillHoles',...
          'Position',[xspot,450,70,25],...
          'Callback',@FillHolesbutton_Callback);
% hSegment = uicontrol('Style','pushbutton','String','Re-segment',...
%           'Position',[xspot,400,70,25],...
%           'Callback',@segmentbutton_Callback);
hSaveImage = uicontrol('Style','pushbutton',...
          'String','SaveImage',...
          'Position',[xspot,300,70,25],...
          'Callback',@savebutton_Callback); 
%    ha = axes(ax,'Units','Pixels','Position',[50,60,200,185]); 
% align([hAdd,hDelete,hSplit,hSaveImage,hpopup,htext,hFillHoles,hRemoveArea],'Center','None');
      
      
      
      
f.Visible = 'on'   ;
If = zeros(2048,2048);
   f.Units = 'normalized';
   htext.Units = 'normalized';
   hpopup.Units = 'normalized';
   hDelete.Units = 'normalized';
   hSplit.Units = 'normalized';
   hSaveImage.Units = 'normalized';
   hAdd.Units = 'normalized';
   hRemoveArea.Units = 'normalized';
   hFillHoles.Units = 'normalized';

% ha = axes('Units','Pixels','Position',[50,60,200,185]); 
% mm = subplot(1,2,1);
mm = axes;
mm.Units ='pixels';
Position = [20 40 990 990];
% Position = [0.05 0.05 0.9 0.9];
mm.Position = Position;
mm.Units = 'normalized';
Position = mm.Position;
himg = imshow(If);


ff= figure;
ff.Visible ='off';
ff.Position =[1220,400,600,600];
ffposition= ff.Position;

set(f,'KeyPressFcn',@keypress);
end


function keypress(fig_obj,~)
global   A If savename
key = get(fig_obj,'CurrentKey');

switch key
    case 'v'
        addareabutton_Callback([],[])
    case 'r' 
        removeareabutton_Callback([],[])
    case 'd'
        deletebutton_Callback([],[])
    case 'f'
        FillHolesbutton_Callback([],[])
    case 'w'
        nextpos([],[])
    case 'q'
        prevpos([],[])
    case 's'
        savebutton_Callback([],[])
    case 'e'
        edgeclear_Callback([],[])
end
end



function nextpos(~,~)
global If DICimage savename pvalue SegLocation HOURS

% Determine the selected data set.

isitthere = strcmp(HOURS,pvalue);
idx = find(isitthere ==1);

    if idx ==1
        idx = idx+1;
    elseif idx == length(HOURS)
    %     idx=idx-1;
    else
        idx=idx+1;
    %     idx=idx-1;
    end
pvalue = HOURS{idx};
setSceneAndUpdate
end

function prevpos(~,~)
global If DICimage savename pvalue SegLocation HOURS
% Determine the selected data set.

isitthere = strcmp(HOURS,pvalue);
idx = find(isitthere ==1);

    if idx ==1
%         idx = idx+1;
    elseif idx == length(HOURS)
        idx=idx-1;
    else
%         idx=idx+1;
        idx=idx-1;
    end
pvalue = HOURS{idx};
setSceneAndUpdate
end


function setSceneAndUpdate
global SegLocation pvalue If DICimage savename MRNAimage IMGlocation
%try for 3G cells first
cd (SegLocation)
SEGfilelist = dir(strcat('*sm*',pvalue,'_*.mat'));
if ~isempty(SEGfilelist)
    load(char(SEGfilelist.name))
    If = segLogical; 
else

SEGfilelist = dir(strcat('*sm*',pvalue,'_*594*.tif'));
if ~isempty(SEGfilelist)
If = imread(char(SEGfilelist.name));
If(If>0)=1;
end
end

%try for C2C12 cells next
if isempty(SEGfilelist)
SEGfilelist = dir(strcat('*C2_',pvalue,'*.mat'));
if ~isempty(SEGfilelist)
    load(char(SEGfilelist.name))
    If = segLogical; 
else
SEGfilelist = dir(strcat('*C2_',pvalue,'*594*.tif'));   
    if ~isempty(SEGfilelist)
        If = imread(char(SEGfilelist.name));
        If(If>0)=1;
    end
end
end

%added on march 2 2016
if isempty(SEGfilelist)
 SEGfilelist = dir(strcat('*sm*',pvalue,'*.mat'))  
 if ~isempty(SEGfilelist)
    load(char(SEGfilelist.name))
    If = segLogical; 
 end
end

DICfilelist = dir(strcat('*DIC*',pvalue,'.tif'));

savename = char(SEGfilelist.name);
DICimage = imread(char(DICfilelist.name));

updateImage
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


function segmentbutton_Callback(source,eventdata) 

global A  strelsize strelsizesub sigma kernelgsize mask_em thresh pvalue If
%%%parameters
strelsize           =    3;      %3
strelsizesub        =    2;     %2
sigma               =    60;    %60     
kernelgsize         =    250;    %250


position = pvalue;
toughsegmentationforFISHembed(position)
% 
% hSegment = uicontrol('Style','pushbutton','String','Re-segment',...
%           'Position',[xspot,400,70,25],...
%           'Callback',@segmentbutton_Callback);

end


 function FillHolesbutton_Callback(source,eventdata) 
 global If AnnotationsDir
 % choose cell
%       [cellx,celly] = ginput(1);
       % construct a polygon to add

    Ig = imfill(If,'holes'); 
      If = Ig;
    updateImage
 end



 function addareabutton_Callback(source,eventdata) 
 global If AnnotationsDir
 % choose cell
%       [cellx,celly] = ginput(1);
       % construct a polygon to add

       [polyx,polyy] = ginput();
       M  = zeros(1,length(polyx)*2);
      M(1:2:end) = polyx;
      M(2:2:end) = polyy;
      zeroImage = zeros(size(If));
      zeroImage = insertShape(zeroImage,'FilledPolygon',M,'LineWidth',6,'Color',[1 1 1]);
      zerogray = rgb2gray(zeroImage);
      If(zerogray>0) = 1;
    updateImage
 end
   
 
  function removeareabutton_Callback(source,eventdata) 
 global If AnnotationsDir
 % choose cell
%       [cellx,celly] = ginput(1);
       % construct a polygon to add
      [polyx,polyy] = ginput();
      

      M  = zeros(1,length(polyx)*2);
      M(1:2:end) = polyx;
      M(2:2:end) = polyy;
      zeroImage = zeros(size(If));
      zeroImage = insertShape(zeroImage,'FilledPolygon',M,'LineWidth',6,'Color',[1 1 1]);
      zerogray = rgb2gray(zeroImage);
      If(zerogray>0) = 0;
    updateImage
   end
 
function deletebutton_Callback(source,eventdata) 
   global stats If AnnotationsDir
   % Display mesh plot of the currently selected data.
     [cellxx,cellyy] = ginput();
      output.cellx=cellxx;
      output.celly=cellyy;
       

      
    for j = 1:length(cellxx)  
        celly = cellyy(j);
        cellx = cellxx(j);
        for i = 1:length(stats)

        m = stats(i).PixelIdxList;
        index = sub2ind(size(If),round(celly),round(cellx));

        if isempty(find(m==index))
        else
        If(m)=0;
        end
        end
    end
    updateImage
end


function edgeclear_Callback(source,eventdata) 
   global stats If AnnotationsDir
   % Display mesh plot of the currently selected data.

       

      
        celly = [ones(1,2048)   1:2048         ones(1,2048).*2048  1:2048];
        cellx = [1:2048         ones(1,2048)    1:2048              ones(1,2048).*2048];
        for i = 1:length(stats)

        m = stats(i).PixelIdxList;
        index = sub2ind(size(If),round(celly),round(cellx));

        if sum(ismember(m,index))<1
        else
        If(m)=0;
        end
        end
 
    updateImage
end


 
function splitbutton_Callback(source,eventdata) 
   global If AnnotationsDir Pvalue
   % Display contour plot of the currently selected data.
      [x,y] = ginput(2);

lines = int32([x(1) y(1) x(2) y(2)]);
If = insertShape(uint8(If),'Line',lines,'LineWidth',6,'Color',[0 0 0]);
If = rgb2gray(If);
updateImage
   end 
 

function popup_menu_Callback(source,eventdata) 
global If DICimage savename pvalue AnnotationsDir SegLocation
% Determine the selected data set.
 str = source.String;
 val = source.Value;
 pvalue = char(str{val});


setSceneAndUpdate
end

function updateImage
global If DICimage savename stats Position mm AnnotationsDir pvalue ff f ffposition
CC = bwconncomp(If);
L = labelmatrix(CC);
stats = regionprops(L,'PixelIdxList');

for i=1:length(stats)
    area= length(stats(i).PixelIdxList);
    m = stats(i).PixelIdxList;
    if area>10000
    else
        If(m) = 0;
    end
end

axes(mm);
children = findobj(mm,'Type','image');
delete(children);
mm.Position = Position;
mm.NextPlot = 'replace';
rgb = label2rgb(L,'lines','k','shuffle');hold on
mm.CLimMode ='auto';
% immy = zeros(size(L));
% immy(L>0) = 255;
% prim = bwperim(immy);
% se = strel('disk',4);
% prim = imdilate(prim,se);
% imgone = rgb(:,:,1);
% imgone(prim) = 255;
% rgb(:,:,1) = imgone;

    hrgb = imagesc(rgb);
    hrgb.AlphaData=0.8;
    titledispname = savename;
    t = regexp(titledispname,'_');
    titledispname(t) = '-';
    title(titledispname);

channelimg = double(DICimage);
cimgline = reshape(channelimg,[1 size(channelimg,1).*size(channelimg,2)]);
prcntl = prctile(cimgline,2);
lprcntl = prctile(cimgline,0.5);
channelimg = uint8(((channelimg-lprcntl)./prcntl).*255);
channelimg(channelimg == 255) =254;
% DICimage = channelimg;
immy = zeros(size(L));
immy(L>0) = 255;
prim = bwperim(immy);
se = strel('disk',2);
prim = imdilate(prim,se);
channelimg(prim) = 255;

j = imagesc(channelimg);
cmap = colormap(gray(255));
cmap(255,:)=[1 0 0];
colormap(cmap);
j.AlphaData = 0.8;
mm.CLim = [0 256];
mm.NextPlot = 'replace';

figure(ff)
ff.Visible = 'on';
if ~isempty(dir(AnnotationsDir))
cd (AnnotationsDir)
end
scene = 'scene*00';
scene(end-1:end) = pvalue(end-1:end);
filelist = dir(strcat('*',scene,'*.tif'));
annfilenames = char({filelist.name});

if ~isempty(annfilenames)
ann = loadUpFinalImageOfStack(annfilenames);
imshow(ann)
else
scene = 's00';
scene(end-1:end) = pvalue(end-1:end);
filelist = dir(strcat('*',scene,'*.jpg'));
annfilenames = char({filelist.name});   
if ~isempty(annfilenames)
ann = imread(annfilenames);
else 
ann = zeros(512,512);
end
imagesc(ann);
h=gca;
end
ff.Position = ffposition;

if isempty(annfilenames)
ann = zeros(size(If));
imshow(ann);hold on
text(500,500,'No Movie Or No Annotations','FontSize',40,'Color',[1 1 1]);hold off
end

imgff = gca;
imgff.YDir = 'reverse';
figure(f)

end


function savebutton_Callback(source,eventdata)
global A If savename
savethatimage(savename,If,A)
end

function savethatimage(savename,If,A)
cd (strcat(A,'autoseg\'));
% imwrite(If,strcat(savename(1:end-4),'.tif'),'Tiff');

segLogical = logical(If);
if isempty(dir(strcat(savename(1:end-4),'.mat')))
save(strcat(savename(1:end-4),'.mat'),'segLogical');
else
save(strcat(savename(1:end-4),'.mat'),'segLogical','-append');
end
disp(strcat('image ... ',savename,'...saved'));  


% segLogical = logical(If);
% if ~isempty(dir(strcat(savename(1:end-4),'.mat')))
% whostruct = whos('-file',strcat(savename(1:end-4),'.mat'));
% if length(whostruct)>1
% save(strcat(savename(1:end-4),'.mat'),'segLogical','-append');
% end
% else
% save(strcat(savename(1:end-4),'.mat'),'segLogical');
% disp(strcat('image ... ',savename,'...saved'));  
% end

end




% Tough Segementation


%% annotations
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

