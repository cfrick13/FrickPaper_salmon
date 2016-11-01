%% smfishflatbkg (dividing of experiment images by flatfield reference images)
% divide
function smfishflatbkgPARALLELnew(Date)
%specify directory
mfile = mfilename('fullpath');
[~,b] = regexp(mfile,'FrickPaperData');
mfiledir = mfile(1:b+1);
parentdir = mfiledir;
A = parentdir;

for BB = Date
B = BB{1};
B = strcat(B,' smFISH');
  
    
C = '\EXPORTS';
H = 'reference';
I = '\FLATFIELD';


REF = strcat(A,B,C);
EXP = REF;
SAV = strcat(A,B,'\',I);
mkdir(SAV)

SAVENAME = 'flat';

        cd (EXP)
        expfilelist = dir(strcat('*z10*.tif'));
        chansinlist = sort(findNumberOfVarsInListP(expfilelist,'(Alexa Fluor 594|Alexa Fluor 647|DIC)'));
        cidx = cellfun(@isempty,chansinlist,'UniformOutput',1);
        chansinlist = chansinlist(~cidx);


for chan = 1:length(chansinlist)
%load flatfield reference image for one channel
    imgback = loadReferenceStack(chansinlist,chan,B,H,REF,SAV);
 
%load the experimental images for one channel
cd (EXP)
filelist = dir(strcat('*',char(chansinlist(chan)),'*.tif'));
cfile = {filelist.name};

    parfor i = 1:length(cfile)
%     for i = 1:length(cfile)
         img = imread(char(cfile{i}));

         dimg = double(img);
         dimgback = double(imgback);

         %new
         backshape = reshape(dimgback,[1 size(dimgback,1).*size(dimgback,2)]);
         backsort = sort(backshape(~isnan(backshape)));
         backtop = backsort(round(length(backsort).*0.999):round(length(backsort).*0.9999));
         dimgbacknorm = dimgback./double(mean(backtop));

         flat = dimg./dimgbacknorm;

         file = strcat(SAVENAME,'_',char(cfile{i}));
         saveThatCorrectedFile(file,SAV,flat,REF)
    end
end
cd('D:\Users\zeiss\Documents\MATLAB\')

end
end


function imgback = loadReferenceStack(chansinlist,chan,B,H,REF,SAV)
    cd (REF);
    filelist = dir(strcat('*',H,'*',char(chansinlist(chan)),'*z10*','*.tif'));
    if isempty(filelist)
        filelist = dir(strcat('*',H,'*',char(chansinlist(chan)),'*.tif'));
    end
    
    cfile = {filelist.name};

%%% determinesize to increase speed
testimg = imread(char(cfile{1}));
imgsize = size(testimg);
clear testimg
%%%

    imgarray = zeros(imgsize(1),imgsize(2),length(cfile));
    for i = 1:length(cfile)
         imgpre = double(imread(char(cfile{i})));

         medpre = nanmedian(nanmedian(imgpre));
         imgpre(imgpre>5*medpre) = NaN;

         %%%
         if sum(sum(isnan(imgpre)))>0
             %%%
             jed = zeros(size(imgpre));
             jed(isnan(imgpre)) = 1;
             jeddil = imdilate(jed,strel('disk',4));
             imgpre(jed==1) = nanmedian(imgpre(jeddil==1));
             %%%
             else
         end
         imgarray(:,:,i) = imgpre;
    end
    
dimgarray = double(imgarray);
med = nanmedian(dimgarray,3);
sig = 30;
kgsize = 60;
imgmed = gaussianBlurz(med,sig,kgsize);%low pass

file = strcat(B,'_',H,'_flat_corr_',char(chansinlist(chan)),'.tif');
saveThatFile(file,SAV,imgmed,REF)
imgback = imgmed;
end

function saveThatFile(file,SAV,med,REF)
    cd (SAV)
    ims = uint16(med);
    imwrite(ims,file,'tiff');
    cd (REF)
end
function saveThatCorrectedFile(file,SAV,flat,REF)
    filename = file;
    disp(filename)
    cd (SAV)
    ims = uint16(flat);
    imwrite(ims,file,'tiff');
    cd (REF)
end
     
function bw = gaussianBlurz(im,sigma,kernelgsize,varargin)

filtersize = [kernelgsize kernelgsize];
kernelg = fspecial('gaussian',filtersize,sigma);

% image filtering
gFrame = imfilter(im,kernelg,'repl');

if ~isempty(varargin)
    bw=gFrame.*uint16(varargin{1}>0);
else
    bw=gFrame;
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