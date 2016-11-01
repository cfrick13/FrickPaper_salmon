function mstackTOdottedWin(Date,director)
%saves a logical matrix of mRNA locations for all scenes in experiment
global   A   zs AnnotationsDir  

ksize =5;
if isempty(director)
    director = strcat('\FISHareaNewTHRESH',num2str(ksize));
else
    director = strcat(director,num2str(ksize));    
end

%%%%%%%%%%%%%%%%%%%%%%%%
% Date = '2015_12_11';
% for Datez = {'2016_02_09','2016_02_11'};
for BB = Date;
    
    Datez = char(BB);
  

%%%%%%%%%%%%%%%%%%%%%%%%
    % close all
    
    mfile = mfilename('fullpath');
    [~,b] = regexp(mfile,'FrickPaperData');
    mfiledir = mfile(1:b+1);
    parentdir = mfiledir;
    A = strcat(parentdir,Datez,' smFISH\');
        
    AnnotationsDir = strcat(A,'\ANNOTATIONS');
    SegLocation = (strcat(A,'\autoseg\'));
    ImgLocation = strcat(A,'\FLATFIELD\');
    MRNALocation = strcat(A,director);
    DottedLocation = strcat(A,'\dotted\');
    mkdir(DottedLocation)
    cd(ImgLocation)



    primarylist = dir('*.tif');
    HOURS = findNumberOfVarsInList(primarylist, 'p[0-9]+');
    channels = findNumberOfVarsInList(primarylist,'(Alexa Fluor 594|Alexa Fluor 647)');
    channn = cellfun(@isempty,channels);
    chann = channels(~channn);
    zs = findNumberOfVarsInList(primarylist,'z[0-9]+');


disp(Datez)
    for chanstr = chann
        fluorchannel = char(chanstr);
%         for i = 1:length(HOURS)
        parfor i = 1:length(HOURS)
            pvalue = char(HOURS{i});
            setSceneAndUpdate(pvalue,fluorchannel,SegLocation,ImgLocation,DottedLocation,MRNALocation) %dotted file is saved within this function
        end
    end
end
end





function setSceneAndUpdate(pvalue,fluorchannel,~,ImgLocation,DottedLocation,MRNALocation)
DotsOn =1;
%try for 3G cells first
cd (ImgLocation)
% MRNAimageStack = imread(char(MRNAfilelist.name));
[details,MRNAimageStack] = makeimagestack(pvalue,fluorchannel);
dotted = loadccmrna(details,pvalue,fluorchannel,DotsOn,MRNAimageStack,ImgLocation,MRNALocation);
cd(DottedLocation)
dottedfilename = strcat('dotted','-',details.species,'-',pvalue,'.mat');

disp(strcat(pvalue,'-',fluorchannel,'-',details.species))
save(dottedfilename,'dotted');

cd(ImgLocation)
% updateImage
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
stack = [];
i=1;
%     for cfile = {imagefilelist.name}
%         stack(:,:,i) = imread(char(cfile));
%         i=i+1;
%     end
    imf = char(imagefilenames{1});
%     [a,b] = regexp(imf,'(594_snail|594_smad7|594_pai1|647_smad7|594_pmepa1|594_tieg|594_bhlhe40|647_snail|647_wnt9a|647_ctgf|594_wnt9a|594_ctgf)');
    [a,b] = regexp(imf,'(594_snail|594_smad7|594_pai1|647_pai1|647_smad7|594_pmepa1|594_tieg|594_bhlhe40|647_snail|647_wnt9a|647_ctgf|594_wnt9a|594_ctgf)');
    [d,e] = regexp(imf,'(off|low|medlow|med|high|0dot00|0dot01|0dot02|0dot03|0dot04|0dot07|2dot40)');
    [f,g] = regexp(imf,'[0-9]hr');
    [c,~] = regexp(imf,'_');
%     imf(c) = '-';
    details.species = imf(a:b);
    details.dose = imf(d:e);
    details.timepoint = imf(f:g);
    
    
end
function dotted = loadccmrna(details,pvalue,~,DotsOn,~,~,MRNALocation)
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
