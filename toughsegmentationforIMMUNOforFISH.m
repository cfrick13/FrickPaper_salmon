function  toughsegmentationforIMMUNOforFISH(Date)

% for DATES = {'2015_01_30 C2C12 immuno timecourse old antibody'}
for DATES = Date

% {'2015_03_06 smFISH'}
Ab = 'D:\Users\zeiss\Pictures\Frick\';
% Ab = 'F:\';
% for DATES = Ba

% for DATES = {'2015_03_06 smFISH'}

B = char(DATES);
B = strcat(B,' smFISH');

A = strcat(Ab,B,'\');
B = 'FLATFIELD\';
C = 'Thresholds7';
mkdir(strcat(A,'autoseg\'));
choosefiles(A,B,C)
end
end

function choosefiles(A,B,Cdir)

cd (strcat(A,B));
primarylist = dir('*Alexa Fluor 594_*');
% primarylist = dir('*_Alexa Fluor 594_*');
HOUR = findNumberOfVarsInListP(primarylist, 'p[0-9]+');
HOURS = sort(HOUR);

for pvaluecell = HOURS
    pvalue = char(pvaluecell);
    disp(pvalue);

cd(strcat(A,Cdir))
file = strcat('*focus*',pvalue,'_*');
filelist = dir(file);
cfile = {filelist.name};
%load focusPoint
load(char(cfile),'focusPoint');
z = num2str(focusPoint);
zvalue = '00';
if length(z)>1
zvalue(end-1:end) = z;
else
zvalue(end) = z;
end

[IfEGFP,savenameEGFP] =  segmentSpecificChannel('Alexa Fluor 594',pvalue,A,B,zvalue);
   


%%
cd (strcat(A,B));
fileDic = dir(strcat('*',pvalue,'-*DIC*z',zvalue,'*.tif'));
if isempty(fileDic)
    stophere=2;
end
DICimage = imread(char(fileDic.name));
filenamezDIC = char(fileDic.name);
savenameDIC = strcat(filenamezDIC(1:18),'_DIC_',pvalue);

% DAPIperim = bwperim(IfDAPI);
EGFPperim = imdilate(bwperim(IfEGFP),strel('disk',2));

C = zeros([size(EGFPperim,1) size(EGFPperim,2) 3]);
% C(:,:,1) = DAPIperim.*4000;
C(:,:,2) = EGFPperim.*4000;

% CC = bwconncomp(IfDAPI);
% L = labelmatrix(CC);
% rgbDAPI = label2rgb(L,'winter','k');
% stats = regionprops(L,'BoundingBox','PixelIdxList');

CC = bwconncomp(IfEGFP);
L = labelmatrix(CC);
rgbEGFP = label2rgb(L,'parula','k','shuffle');
stats = regionprops(L,'BoundingBox','PixelIdxList');

figure(10)
hh = imshow(rgbEGFP);hold on
hh.AlphaData =1;
% h = imshow(rgbDAPI);hold on
% h.AlphaData=0.5;
j = imagesc(DICimage);hold on
j.AlphaData = 0.8;
colormap('gray')
meanDIC = mean(mean(DICimage));
lowDIC = round(meanDIC./2);
highDIC = round(meanDIC.*2);
hax = gca;
hax.CLim = [lowDIC highDIC];
h = imagesc(C);
h.AlphaData=0.3;
hold off
drawnow

stophere=1;

% savethatimage(savenameDAPI,uint8(IfDAPI),A)
savethatimage(savenameEGFP,uint8(IfEGFP),A)
savethatimage(savenameDIC,uint16(DICimage),A)

end
end

function [If,savename] =  segmentSpecificChannel(CHANNEL,pvalue,A,B,zvalue)
 
% DAPI segmented image
% cd (strcat(A,B));
% file = dir(strcat('*','*',pvalue,'*',CHANNEL,'*z',zvalue,'*'));
% FocusImage = imread(char(file.name));
% filenamez = char(file.name);
% savename = strcat(filenamez(1:18),'_',pvalue,'_*',CHANNEL);
% disp(savename);



cd (strcat(A,B));
file = dir(strcat('*','*',pvalue,'-*',CHANNEL,'*z','*'));
zlength = length(file);
FocusImages = zeros(2048,2048,zlength);

for z =1:zlength
    zvalue ='00';
   znum = num2str(z);
   if length(znum)>1
       zvalue(end-1:end) = znum;
   else
       zvalue(end) = znum;
   end
   file = dir(strcat('*','*',pvalue,'-*',CHANNEL,'*z',zvalue,'*'));
   FocusImages(:,:,z) = imread(char(file.name));
end

FocusImage = max(FocusImages,[],3); %maximum intensity projection!
filenamez = char(file.name);
savename = strcat(filenamez(1:18),'_',pvalue,'_',CHANNEL);
disp(savename);




if strcmp(CHANNEL,'DAPI')
% If = segmentationDAPI(FocusImage);
% If = segmentationDAPI(FocusImage);
elseif strcmp(CHANNEL,'Alexa Fluor 594')
If = segmentationEGFP(FocusImage);
end

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
    if ~isempty(hours)
    HOURS{jjj} = hours;
    jjj=jjj+1;
    end
else
end
end
end
end


function FinalImage=loadStack(FileTif)
% [a,b] = uigetfile;
% FileTif = a;
% cd (b)
InfoImage=imfinfo(FileTif);
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
NumberImages=length(InfoImage);
FinalImage=zeros(nImage,mImage,NumberImages,'uint16');
 
TifLink = Tiff(FileTif, 'r');
for i=1:NumberImages
   TifLink.setDirectory(i);
   FinalImage(:,:,i)=TifLink.read();
end
TifLink.close();
end


function bw = gaussianBlurz(im,sigma,kernelgsize,varargin)

filtersize = [kernelgsize kernelgsize];
kernelg = fspecial('gaussian',filtersize,sigma);

%% image filtering
gFrame = imfilter(im,kernelg,'repl');

if ~isempty(varargin)
    bw=gFrame.*uint16(varargin{1}>0);
else
    bw=gFrame;
end
end

function If = segmentationDAPI(img)
fig=1;
%% parameters
strelsize           =    3;      %3
strelsizesub        =    2;     %2
sigma               =    5;    %60     
kernelgsize         =    10;    %200
mask_em         = zeros(size(img));
% left = 0.009;
% slopedown = 0.007;
% % left = 0.005;
% % slopedown = 0.0045;
left = 0.012;
slopedown = 0.0001;


%% start
imgorig = img;
%smooth initial image
se = strel('disk',strelsize);

%%%
Ie = imerode(img,se);
%%%
sig = 20;
kgsize=20;
ILP = gaussianBlurz(img,sig,kgsize);%low pass
Ie = imerode(ILP,se);
imgt=ILP;
%%%
img=imgt;

Iobr = imreconstruct(Ie,img);
Iobrd = imdilate(Iobr,se);
Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
img = Iobrcbr;
img=imgorig;



% 
%     figure(fig)
%     fig = fig+1;
%     imagesc(img);
%     title('smoothed input image');
    
    
% 
% 
%             %ONLY NECESSARY IF TRYING TO FIND THE PEAKS IN THE SMOOTHED GAUS 
%             
    gaus = gaussianBlurz(imgorig,sigma,kernelgsize);
%            
%     figure(fig)
%     fig = fig+1;
%     hold off
%     bar(bincenters,fraction)
%     ylim([0 0.01])
%     ylim([0 0.003])
%     xlim([0 6000])
%     hold on
%     stem(segthresh,1,'g');
%     hold off


% findpeaksgaus(gaus) % homemade function to find peaks after applying gaussian blur 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sub = double(img) -double(gaus);%%%%%%% key step!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
div = double(gaus)./double(img);
cropz =10;
maxlocale = max(max(div(cropz:end-cropz,cropz:end-cropz)));
d = find(div == maxlocale);
div(d);
positiond = d(1);
[e,f] = ind2sub(size(img),positiond);


minlocale = min(min(sub(cropz:end-cropz,cropz:end-cropz)));
a = find(sub == minlocale);
sub(a);

position = a(1);
[b,c] = ind2sub(size(img),position);

%     figure(fig)
%     fig = fig+1;
%     imagesc(sub);
%     title('subtracted');
%     hold on
%     plot(c,b,'y+');


displaynormalizationfactor(e,f)
%normalizationfactors = double(gaus(b,c))./double(img(b,c));
normalizationfactors = double(gaus(e,f))./double(img(e,f));
normalizationfactor = mean(mean(normalizationfactors));

% figure(fig)
% fig = fig+1;
% imagesc(subgauss);hold on
% plot(c,b,'y+')

scaledgaus = double(gaus)./normalizationfactor;
sub_scale_corr = double(img)-double(scaledgaus);

 
%     figure(fig)
%     fig = fig+1;
%     imagesc(sub_scale_corr);
%     title('subtracted_scaled');


subtractionref = sub_scale_corr;
vec = reshape(subtractionref,size(subtractionref,1)^2,1);
% [numbers,bincenters] = hist(double(vec),0:0.5:500);
% 
% 
% 
% numbers = medfilt1(numbers, 10); %smooths curve
% fraction = numbers./sum(numbers);
% 
% 
% % leftedge = find(fraction > max(fraction.*0.80),1,'first');
% % insideslopedown = find(fraction(leftedge:end) <max(fraction.*0.3),1,'first');
% 
% 
% 
% leftedge = find(fraction > left,1,'first');
% insideslopedown = find(fraction(leftedge:end) < slopedown,1,'first');
% 
% insideslopeup = find(fraction(leftedge+insideslopedown:end) >0.0012,1,'first');
% trough = min(fraction(leftedge+insideslopedown:leftedge+insideslopedown+insideslopeup));
% troughindex = find(fraction(leftedge+insideslopedown:leftedge+insideslopedown+insideslopeup) == trough);
% troughindexrounded = round(median(troughindex));
% 
% % threshlocation = bincenters(leftedge+insideslopedown+troughindexrounded);
% 
% threshlocation = bincenters(leftedge+insideslopedown);
% % 
% %     figure(fig)
% %     fig=fig+1;
% %     bar(bincenters,fraction);hold on
% %     xlim([5 150])
% %     ylim([0 0.01])
% %     stem ([threshlocation threshlocation],[0 1]);hold off
% 
% subtractionthreshold = threshlocation;
% 
% if size(subtractionthreshold,1)==size(subtractionthreshold,2)
% else
%      subtractionthreshold = mean(threshlocation);
% end
% 

sortedvec = sort(vec);
toppest=0.999;
topper=0.97;
vecTOP = nanmean(sortedvec(round(length(sortedvec).*topper):round(length(sortedvec).*toppest)));

subtractionthreshold = vecTOP;

subtracted = sub_scale_corr-subtractionthreshold;
subzero = (subtracted<0);
subtractedzero = subtracted.*(~subzero);

%     figure(fig)
%     fig = fig+1;
%     imagesc(subtractedzero,[0 20]);
%     title('subtracted_scaled_zeroed');
% 
% %     a = find(subtractedzero<4);
% % submax = zeros(size(subtractedzero));
% % subtractedzero(a)=0;
% 

% se = strel('disk',strelsizesub);
% Ie = imerode(If,se);

% If = subtractedzero;
% Ie = imerode(If,se);


Ie = subtractedzero;
a = find(Ie>0);
submax = zeros(size(Ie));
Ie(a)=50;

se = strel('disk',2);
If = imclose(Ie,se);
% se = strel('disk',1);
% If = imerode(Ig,se);
% 

% imgt = -double(imgt);
% imgt(~(If==50)) = -Inf;


sigma               =    2;    %60     
kernelgsize         =    40;    %200
gaus = gaussianBlurz(gaus,sigma,kernelgsize);
sigma               =    2;    %60     
kernelgsize         =    40;    %200
gaus = gaussianBlurz(gaus,sigma,kernelgsize);
sigma               =    2;    %60     
kernelgsize         =    40;    %200
gaus = gaussianBlurz(gaus,sigma,kernelgsize);
sigma               =    2;    %60     
kernelgsize         =    40;    %200
gaus = gaussianBlurz(gaus,sigma,kernelgsize);
sigma               =    2;    %60     
kernelgsize         =    40;    %200
gaus = gaussianBlurz(gaus,sigma,kernelgsize);
imgt = -double(gaus);
imgt(~(If>0)) = -Inf;

L=watershed(imgt);
If = L>1;


stophere=1;


end


function savethatimage(savename,If,A)
cd (strcat(A,'autoseg\'));
imwrite(If,strcat(savename,'.tif'),'Tiff');
cd ..
end

function displaynormalizationfactor(b,c)
disp(strcat(num2str(b),'_',num2str(c)));
end

function If = segmentationEGFP(img)
fig=1;
%% parameters
left = 0.004;
slopedown = 0.003;

firststrel = 30;
sigmafirst = firststrel.*3;
kernelgsizefirst = firststrel.*6;

%% start
imgorig = img;

img = wiener2(img,[5 5]);
se =strel('disk',2);
Ie = imerode(imgorig,se);
Iobr = imreconstruct(Ie,img);
Iobrd = imdilate(Iobr,se);
Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
IobrcbrF = imcomplement(Iobrcbr);
gaus = double(IobrcbrF);

se =strel('disk',firststrel);
Ie = imerode(gaus,se);
Iobr = imreconstruct(Ie,gaus);
Iobrd = imdilate(Iobr,se);
Iobrcbr = imreconstruct(imcomplement(Iobrd),imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
gaus = double(Iobrcbr);


sigma = sigmafirst;
kernelgsize = kernelgsizefirst;
gaustwo = gaussianBlurz(double(gaus),sigma,kernelgsize);

sub = double(gaus) -double(gaustwo);%%%%%%% key step!
b = find(sub == min(min(sub)));
rattio = gaustwo(b)./gaus(b);
gaustwocorr = gaustwo./rattio;
sub_scale_corr = double(gaus) - double(gaustwocorr);




subtractionref = sub_scale_corr;
vec = reshape(subtractionref,size(subtractionref,1)^2,1);
[numbers,bincenters] = hist(double(vec),0:0.5:10000);



numbers = medfilt1(numbers, 10); %smooths curve
fraction = numbers./sum(numbers);

mf = max(fraction);

%%%%%%%%%%%%%%%%%%%%
left=0.3*mf;
slopedown=0.2*mf;
%%%%%%%%%%%%%%%%%%%%%

leftedge = find(fraction > left,1,'first');
insideslopedown = find(fraction(leftedge:end) < slopedown,1,'first');
insideslopeup = find(fraction(leftedge+insideslopedown:end) >0.0012,1,'first');
trough = min(fraction(leftedge+insideslopedown:leftedge+insideslopedown+insideslopeup));
troughindex = find(fraction(leftedge+insideslopedown:leftedge+insideslopedown+insideslopeup) == trough);
troughindexrounded = round(median(troughindex));



threshlocation = bincenters(leftedge+insideslopedown);
% 
% figure(22)
%     bar(bincenters,fraction);hold on
%     xlim([-500 1000])
%     ylim([0 0.1])
% stem ([threshlocation threshlocation],[0 1]);hold off
% drawnow

subtractionthreshold = threshlocation;

if size(subtractionthreshold,1)==size(subtractionthreshold,2)
else
     subtractionthreshold = mean(threshlocation);
end

% subtractionthreshold = graythresh(subtractionref);

subtracted = sub_scale_corr-subtractionthreshold;
subzero = (subtracted<0);
subtractedzero = subtracted.*(~subzero);


Ie = subtractedzero;
a = find(Ie>0);
submax = zeros(size(Ie));
Ie(a)=50;

Ih = Ie>0;
Igclose = imclose(Ih,strel('disk',30));
Igclosemax = imclose(Ih,strel('disk',80));
Igcopenmax = imopen(Igclosemax,strel('disk',10));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
se = strel('disk',12);
Igcopenmax = imdilate(Igcopenmax,se);

Igcopen = imopen(Igclose,strel('disk',2));
Igcofill = imfill(Igcopen,'holes');
Igcfopen = bwareaopen(Igcofill,5000);
Igcfopendil = imdilate(Igcfopen,strel('disk',5));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
se = strel('disk',12);
Igcfopendil = imdilate(Igcfopendil,se);

sigma=20;
kernelgsize=40;
gaus = gaussianBlurz(IobrcbrF,sigma,kernelgsize);

sigma=40;
kernelgsize=80;
gaus = gaussianBlurz(gaus,sigma,kernelgsize);

imgt = -double(gaus);
% imgt(~(Igcfopen>0)) = -Inf;
imgt(~(Igcopenmax>0)) = -Inf;

L=watershed(imgt);

L(Igcfopendil<1) = 0;
imagesc(L)
colormap parula
If = L>1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% se = strel('disk',8);
% If = imdilate(If,se);



% If = imfill(Ih,'holes');
%  imagesc(If);
%  
%     figure(fig)
%     fig=fig+1;
%     imagesc(If);
%     title('subtractedzeromax');


%     
% gausmin = gaussianBlurz(img,40,80);
% gausmin = -single(gausmin);
% 
% gausmin(~logical(If)) = -Inf;
% Lgaus = watershed(gausmin);
% rgbgaus = label2rgb(Lgaus,'jet' ,[.6 .6 .6]);

% se = strel('disk',3);
% Iwater = imerode(If,se);
% D = bwdist(~Iwater,'euclidean');
% Dgaus  = gaussianBlurz(D,40,40);
% D=-Dgaus;
% D(~If) = -Inf;
% L = watershed(D);
% rgb = label2rgb(L,'jet' ,[.5 .5 .5]);

%     figure(fig)
%     fig=fig+1;
%         fig =1;
% imshow(rgbgaus)
stophere=1;


% Lgausplus = Lgaus+1;
% Lgausplus(Lgausplus==2)=1;
% Lgausplus(Lgausplus==1)=0;
% Lgausplus(Lgausplus>0)=1;

% time = settimecharacter(frames);

% imwrite(uint8(Ie),strcat(scenename,'_','t',time,'_c4_flat.tif'));
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