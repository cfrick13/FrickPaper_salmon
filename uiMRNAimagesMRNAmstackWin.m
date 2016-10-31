function uiMRNAimagesMRNAmstackWin
global CLOW CHIGH DotsOn hchannels hslice hpopup If DICimage hproject hdots stats MRNALocation project A Position mm zs AnnotationsDir pvalue ff f SegLocation HOURS ImgLocation fluorchannel zslice

close all
%%%%%%%%%%%%%%%%%%%%%%%%
% Date = '2015_12_11';
Date = '2016_02_09';
%%%%%%%%%%%%%%%%%%%%%%%%
CLOW = 400;
CHIGH = 1600;

DotsOn =1;
stats=struct();
DICimage = zeros(2048,2048);
% close all
% A = strcat('D:\Users\zeiss\Pictures\Frick\',Date,' smFISH\');
% A = strcat('D:\Users\zeiss\Pictures\Frick\',Date,' smFISH');
A = strcat('F:\FrickPaper\',Date,' smFISH\');
% A = strcat('\Users\frick\Documents\Goentoro_Lab\DATA\current\',Date,' smFISH\');
% AnnotationsDir = strcat(A,'ANNOTATIONS\');
AnnotationsDir = strcat(A,'ANNOTATIONS');
SegLocation = (strcat(A,'autoseg\'));
ImgLocation = strcat(A,'FLATFIELD\');
% MRNALocation = strcat(A,'FISHareaNewMEDIAN5\');
MRNALocation = strcat(A,'FISHareaNewTHRESH5\');
cd(ImgLocation)



primarylist = dir('*.tif');
HOURS = findNumberOfVarsInList(primarylist, 'p[0-9]+');
channels = findNumberOfVarsInList(primarylist,'(DIC|Alexa Fluor 594|Alexa Fluor 647)');
zs = findNumberOfVarsInList(primarylist,'z[0-9]+');


fluorchannel = channels{1};
zslice = 1;
project = 'slices';


f = figure;
f.Visible ='off';
f.Position =[10,50,1300,1050];
xspot = 1100;

htext = uicontrol('Style','text','String','Choose Position to Correct',...
          'Position',[xspot,550,60,15]);
hpopup = uicontrol('Style','popupmenu',...
          'String',HOURS',...
          'Position',[xspot,500,100,25],...
          'Callback',@popup_menu_Callback);
      
hchannels = uicontrol('Style','popupmenu',...
          'String',channels',...
          'Position',[xspot,600,100,25],...
          'Callback',@channels_Callback);
      
hslice = uicontrol('Style','popupmenu',...
          'String',zs',...
          'Position',[xspot,700,100,25],...
          'Callback',@slice_Callback);
      
hproject = uicontrol('Style','popupmenu',...
          'String',{'slices','projection'},...
          'Position',[xspot,400,100,25],...
          'Callback',@project_Callback);
      
hdots = uicontrol('Style','popupmenu',...
          'String',{'Dots ON','Dots OFF'},...
          'Position',[xspot,300,100,25],...
          'Callback',@dots_Callback);
 
      
f.Visible = 'on'   ;
If = zeros(2048,2048);
   f.Units = 'normalized';
   htext.Units = 'normalized';
   hpopup.Units = 'normalized';
   hchannels.Units = 'normalized';
   hDelete.Units = 'normalized';
   hSplit.Units = 'normalized';
   hSaveImage.Units = 'normalized';
   hAdd.Units = 'normalized';
   hRemoveArea.Units = 'normalized';
   hFillHoles.Units = 'normalized';
   hslice.Units = 'normalized';
   hproject.Units = 'normalized';

% ha = axes('Units','Pixels','Position',[50,60,200,185]); 
% mm = subplot(1,2,1);
mm = axes;
mm.Units = 'pixels';
Position = [20 20 1000 1000];
mm.Position = Position;
mm.Units = 'normalized';
Position = mm.Position;
% Position(3:4)= [0.9 0.9];
% Position = [0.05 0.05 0.9 0.9];
% mm.Position = Position;
himg = imagesc(If);
% mm.Position = Position;


ff= figure;
ff.Visible ='off';
ff.Position =[1100,10,1000,1000];

set(f,'KeyPressFcn',@keypress);
end


function keypress(fig_obj,~)
global   CLOW CHIGH A If savename hslice hpopup hproject hchannels hdots
key = get(fig_obj,'CurrentKey');

switch key
    case 'e'
        hproject.Value = 1;
        project_Callback(hproject,[])
    case 'r' 
        hproject.Value = 2;
        project_Callback(hproject,[])
    case 'f'
%         FillHolesbutton_Callback([],[])
    case 'rightarrow'
        nextpos(hpopup,[])
    case 'leftarrow'
        prevpos(hpopup,[])
    case 'uparrow'
        nextslice(hslice,[])
    case 'downarrow'
        prevslice(hslice,[])
    case 's'
        hchannels.Value = 1;
        channels_Callback(hchannels,[])
    case 'd'
        hchannels.Value = 3;
        channels_Callback(hchannels,[])
    case 'i'
        hdots.Value = 1;
        dots_Callback(hdots,[])
    case 'o'
        hdots.Value = 2;
        dots_Callback(hdots,[])
    case 'c'
        prompt = {'clow','chigh'};
        dlg_title = 'set the contrast';
        ctrast = (inputdlg(prompt,dlg_title));
        CLOW = str2num(ctrast{1});
        CHIGH = str2num(ctrast{2});
        updateImage
end
end



function nextpos(source,~)
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
source.Value = idx;
setSceneAndUpdate
end
function prevpos(source,~)
global If DICimage savename pvalue SegLocation HOURS MRNAimage 
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
source.Value = idx;
setSceneAndUpdate
end


function nextslice(source,~)
global If DICimage savename pvalue SegLocation zs zslice

% Determine the selected data set.
zc = num2str(zslice);
zchar = 'z00';
zchar(end-length(zc)+1:end) = zc;

isitthere = strcmp(zs,zchar);
idx = find(isitthere ==1);

    if idx ==1
        idx = idx+1;
    elseif idx == length(zs)
    %     idx=idx-1;
    else
        idx=idx+1;
    %     idx=idx-1;
    end
zslice = zs{idx};
zchar = zs{idx};
zc = str2num(zchar(end-1:end));
zslice=zc;
source.Value = zc;
updateImage
end
function prevslice(source,~)
global If DICimage savename pvalue SegLocation zs zslice

% Determine the selected data set.
zc = num2str(zslice);
zchar = 'z00';
zchar(end-length(zc)+1:end) = zc;

isitthere = strcmp(zs,zchar);
idx = find(isitthere ==1);

    if idx ==1
%         idx = idx+1;
    elseif idx == length(zs)
        idx=idx-1;
    else
%         idx=idx+1;
        idx=idx-1;
    end
zchar = zs{idx};
zc = str2num(zchar(end-1:end));
zslice=zc;
source.Value = zc;
updateImage
end

function project_Callback(source,eventdata)
global zslice project
% Determine the selected data set.
 str = source.String;
 val = source.Value;
project = char(str{val});

updateImage
end

function dots_Callback(source,~)
global zslice project DotsOn dotted
% Determine the selected data set.
%  str = source.String;
 val = source.Value;
DotsOn = val;

if DotsOn ==2
    dotz = dotted;
    dotted = zeros(size(dotted));
else
    setSceneAndUpdate
end
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

function channels_Callback(source,~) 
global fluorchannel
% Determine the selected data set.
 str = source.String;
 val = source.Value;
 fluorchannel = char(str{val});


updateImage
end

function slice_Callback(source,eventdata)
global zslice
% Determine the selected data set.
 str = source.String;
 val = source.Value;
 zchar = char(str{val});
zslice = str2num(zchar(2:end));

updateImage
end


function setSceneAndUpdate
global SegLocation pvalue cellIDs If DotsOn DICimage savename details cellLocations ImgLocation fluorchannel MRNAimageStack CCmrna dotted mRNAinCell centroid
%try for 3G cells first
cd (SegLocation)
cellLocations = [];
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



cd (SegLocation)
DICfilelist = dir(strcat('*DIC*',pvalue,'.tif'));
if isempty(DICfilelist)
   DICimage = zeros(2048,2048); 
else
DICimage = imread(char(DICfilelist.name));
end
cd (ImgLocation)
% MRNAimageStack = imread(char(MRNAfilelist.name));
[details,MRNAimageStack] = makeimagestack(pvalue,fluorchannel);
savename = char(SEGfilelist.name);


dotted = loadccmrna(details,pvalue,fluorchannel,DotsOn,MRNAimageStack);

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
HOURS = sort(HOURS);
end

function updateImage
global CLOW CHIGH cellIDs cellLocations project dotted details mRNAinCell centroid If DICimage savename stats Position mm AnnotationsDir pvalue ff f MRNAimageStack fluorchannel zslice
 
mstack = MRNAimageStack;
z = zslice;
dottedz = dotted;
mRNAinCell;
centroid;




imski = bwperim(If);
imski = imdilate(imski,strel('square',4));
CC = bwconncomp(imski);
L = labelmatrix(CC);

%clear the axes from the previous frame
axes(mm);
children = findobj(mm,'Type','image');
delete(children);
    mm.Position = Position;
    mm.NextPlot = 'replace';



%tile
titledispname = savename;
t = regexp(titledispname,'_');
titledispname(t) = '-';
titlearray{1} = strcat(details.species,'.   .',details.timepoint,'.    .',details.dose);
titlearray{2} = pvalue;
titlearray{3} = num2str(zslice);
title(titlearray,'FontSize',14);

if strcmp(project,'projection')
   mstack = max(MRNAimageStack,[],3); 
   dottedz = sum(dotted,3);
   z=1;
end


%display DIC image to show
if strcmp(fluorchannel,'DIC')

        %make color image
    rgb = label2rgb(L,'jet','k','shuffle');hold on
    hrgb = imshow(rgb);
    hrgb.AlphaData=0.8;


        j = imagesc(DICimage);
        j.AlphaData = 0.8;
        colormap('gray')
        h = gca;
        h.CLim=([median(median(DICimage)).\2 median(median(DICimage)).*2]);
        mm.NextPlot = 'replace';
elseif strcmp(fluorchannel,'Alexa Fluor 594')
    
        
        mrnaimage = mstack(:,:,z);
        dottedimage = logical(dottedz(:,:,z));

    img = zeros(size(L));
    img(dottedimage) = 255;    
        %make color image
    rgb = double(label2rgb(L,'jet','k','shuffle'));hold on
    rgb(:,:,1) = rgb(:,:,1)+img;
    hrgb = imshow(rgb);
    hrgb.AlphaData=0.8;

%         cmap = colormap(gray(255));%65535
%         cmap(end+1,:) = [1 0 0];
%         mrnaimage(dottedimage) = 256;
        j = imagesc(mrnaimage);
        j.AlphaData = 0.8;
        colormap('gray')
        h = gca;
        h.CLim=([CLOW CHIGH]);
%         h.CLim=([200 800]);
        mm.NextPlot = 'replace';    
elseif strcmp(fluorchannel,'Alexa Fluor 647')
%     z = zslice;
%     mrnaimage = MRNAimageStack(:,:,z);
%     j = imagesc(mrnaimage);
%     j.AlphaData = 0.8;
%     colormap('gray')
%     h = gca;
%     h.CLim=([median(median(mrnaimage)).\2 median(median(mrnaimage)).*2]);
%     mm.NextPlot = 'replace'; 
        mrnaimage = mstack(:,:,z);
        dottedimage = logical(dottedz(:,:,z));

    img = zeros(size(L));
    img(dottedimage) = 255;    
        %make color image
    rgb = double(label2rgb(L,'jet','k','shuffle'));hold on
    rgb(:,:,1) = rgb(:,:,1)+img;
    hrgb = imshow(rgb);
    hrgb.AlphaData=0.8;

%         cmap = colormap(gray(255));%65535
%         cmap(end+1,:) = [1 0 0];
%         mrnaimage(dottedimage) = 256;
        j = imagesc(mrnaimage);
        j.AlphaData = 0.8;
        colormap('gray')
        h = gca;
        h.CLim=([CLOW CHIGH]);
%         h.CLim=([200 800]);
        mm.NextPlot = 'replace';    
else
    
    
    %make color image
rgb = label2rgb(L,'jet','k','shuffle');hold on
hrgb = imshow(rgb);
hrgb.AlphaData=0.8;


    j = imagesc(DICimage);
    j.AlphaData = 0.8;
    colormap('gray')
    h = gca;
    h.CLim=([median(median(DICimage)).\2 median(median(DICimage)).*2]);
    mm.NextPlot = 'replace';
end

for i = 1:length(mRNAinCell)
    hold on;
    xy = centroid{i};
text(xy(1),xy(2),num2str(mRNAinCell(i)),'FontSize',14,'Color','r');
end

if ~isempty(cellLocations)
    for i=1:length(cellLocations)
xy = cellLocations{i};
text(xy(1),xy(2),strcat('cell#',cellIDs{i}),'FontSize',8,'Color','y');
    end
end

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




%% Tough Segementation


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





%% image stack loader
function [details,stack] = makeimagestack(pvalue,chan)
imagefilelist = dir(strcat('*',pvalue,'-*',chan,'*.tif'));
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
    imf = char(imagefilenames{1});
    [a,b] = regexp(imf,'(594_snail|594_smad7|594_pai1|594_pmepa1|594_tieg|594_bhlhe40|647_snail|647_wnt9a|647_ctgf|594_wnt9a|594_ctgf)');
    [d,e] = regexp(imf,'(off|low|medlow|med|high|0dot00|0dot01|0dot02|0dot03|0dot04|0dot07|2dot40)');
    [f,g] = regexp(imf,'[0-9]hr');
    [c,~] = regexp(imf,'_');
%     imf(c) = '-';
    details.species = imf(a:b);
    details.dose = imf(d:e);
    details.timepoint = imf(f:g);
    
    
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
% image filtering
lapFrame = imfilter(im,kernel,'repl');
if ~isempty(varargin)
    bw=lapFrame.*uint16(varargin{1}>0);
else
    bw=lapFrame;
end
end

function dotted = loadccmrna(details,pvalue,fluorchannel,DotsOn,stack)
global  ImgLocation MRNALocation
% cd (ImgLocation)
cd (MRNALocation)

channy = details.species;
filelist = dir(strcat('mstack*',channy,'*',pvalue,'*.mat'));
filename = char({filelist.name});
A = load(filename,'mstack');
dottedz=A.mstack;




CC = bwconncomp(dottedz,6);
    px = CC.PixelIdxList;
    pxl = arrayfun(@(x) length(x{1}),px,'UniformOutput',1);
    pxz = pxl>1;
    pxlist = px(pxz);
        CC.NumObjects = length(pxlist);
        CC.PixelIdxList = pxlist;


arealimit=1;
dotted = false(CC.ImageSize);%initialize the variable
        zul = zeros(1,length(CC.PixelIdxList));
        for q = 1:length(CC.PixelIdxList)
            if length(CC.PixelIdxList{q}) < arealimit
            else
            dot = CC.PixelIdxList{q};

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% Make the DOTTED stack that circles mrna for count validation
            [x,y,z] = ind2sub([2048 2048],dot);
            
            meanx = round(nanmean(x));
            meany = round(nanmean(y));
            meanz = round(nanmean(z));
            zul(q) = length(unique(z));
            
            
            if DotsOn ==1;
                if meanx>2045 || meanx <4 || meany<4 || meany>2046
                    dotted(meanx:meanx,meany:meany,meanz)=1;
%                     disp('edge')
                else
                    if zul(q)<6
                        dotted(meanx:meanx,meany:meany,meanz)=1;
                    elseif zul(q)>19
                        dotted(meanx-2:meanx+2,meany-2:meany+2,unique(z))=1;
                    else
                        zq = unique(z);
%                         dotted(meanx-3:meanx+3,meany:meany,1)=1; %!
%                         dotted(meanx:meanx,meany-3:meany+3,1)=1; %1
                        dotted(meanx-2:meanx+2,meany:meany,zq(round(zul(q)./3)))=1; %1/3
                        dotted(meanx-2:meanx+2,meany:meany,zq(round((zul(q).*2)./3)))=1; %2/3
                    end
                end
                

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            end
            

            end
        end
        
            if DotsOn ==1
            SE = strel('disk',4);
            doots = imdilate(dotted,SE);
            dotted = bwperim(doots,8);
            end
       

    
stopehre=1;

end
