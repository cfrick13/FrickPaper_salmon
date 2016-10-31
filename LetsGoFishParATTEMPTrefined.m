function LetsGoFishParATTEMPTrefined(Date,director)

A = 'D:\Users\zeiss\Pictures\Frick\';
% A = Ab;
% for BB = Ba
% for BB = {'2015_01_15 smFISH'};
% for BB = {'2015_01_31 smFISH','2015_03_06 smFISH'};
% for BB = {'2015_01_15 smFISH','2015_01_19 smFISH','2015_01_29 smFISH','2015_01_30 smFISH','2015_01_31 smFISH','2015_03_06 smFISH','2015_03_25 smFISH','2015_03_31 smFISH','2015_04_01 smFISH',}
% for BB = {'2015_01_15 smFISH','2015_01_19 smFISH','2015_01_29 smFISH','2015_01_31 smFISH','2015_03_06 smFISH','2015_03_25 smFISH','2015_03_31 smFISH','2015_04_01 smFISH','2015_05_14 smFISH','2015_07_02 smFISH','2015_07_10 smFISH','2015_08_31 smFISH','2015_09_03 smFISH'};
for BB = Date;
    
    B = BB{1};
    B = strcat(B,' smFISH');
    D = '\FLATFIELD';
    F = director;
    
    k647 = 5;
    k594 = 5;

    for ksize = 5;


        %set directories
        EXP = strcat(A,B,D); 
        SAV = strcat(A,B,F,num2str(ksize(1)));
        mkdir(SAV)
        cd (EXP)

        %determine values present
        filelist = dir(strcat('*','*.tif'));
        
        
        speciestokens = '(594_snail|594_smad7|594_pai1|647_smad7|594_pmepa1|594_pai1|594_tieg|594_bhlhe40|647_wnt9a|647_ctgf|594_wnt9a|594_ctgf|647_smad7|647_tieg|647_bhlhe40|647_snail)';
        CHANSS = findNumberOfVarsInList(filelist,'(Fluor 594|Fluor 647)');
        CHANS = cellfun(@(x) x(end-2:end),CHANSS,'UniformOutput',0);

        %determine number of positions to analyze
        cd (EXP)
        expfilelist = dir(strcat('*DIC*z10*.tif'));
        PVALUES = sort(findNumberOfVarsInListP(expfilelist,'p[0-9]+'));



        disp(strcat('process',B))
        %compute the LoG convolved and thresholded images
            for chans = CHANS
                chan = char(chans);

                    if strcmp(chan,'647')
                    ksize = k647;
                    else
                    ksize = k594;
                    end

%                 for cfile = 1:length(PVALUES)
                parfor cfile = 1:length(PVALUES)
                pvalue = char(PVALUES{cfile});
                cd(EXP)
                [details,fishdotstack] = makeimagestack(EXP,pvalue,chan,ksize);
            [onlysegmented,thresh] = findthresh(fishdotstack,details,SAV);
                CCmrna = bwconncomp(onlysegmented,6);%%%????? is 6 or 18 or 26 better?
                savemrna(CCmrna,details,thresh,SAV)
                end
            end
    end
end
end



function [dottedz,thresh] = findthresh(fishdotstack,details,SAV)
threshmat = zeros(1,size(fishdotstack,3));
for ii = 1:size(fishdotstack,3)
    ff = fishdotstack(:,:,ii);
        [~,cm] =FastPeakFind(ff);
        cmlog = logical(cm);
        cml = cmlog;
    fpvec = ff(cml);
    fvec = ff(~isnan(ff));
    baseliner = median(fvec);
        [Nn,Eedges] = histcounts(fpvec,baseliner:100:10000,'Normalization','probability');
        [NN,~] = histcounts(fvec(fvec>baseliner),baseliner:100:10000,'Normalization','probability');
%         [Nn,Eedges] = histcounts(fpvec,0:100:10000,'Normalization','probability');
%         [NN,~] = histcounts(fvec,0:100:10000,'Normalization','probability'); %orignal
        Nm = Nn+NN;

    roughhistsm = smooth(Nm);
    roucghhistsmo = smooth(roughhistsm);
    roughhistsmoo = smooth(roucghhistsmo);
    [~,locs] = findpeaks(1-roughhistsmoo);
    
    
    histgrad = gradient(roughhistsmoo);
    [~,locz] = findpeaks(histgrad);
%     %figure for checking
%     figure(77)
%     subplot(2,2,1);
%     plot(roughhistsmoo./max(roughhistsmoo));hold on
%     plot(histgrad./max(histgrad))
%     ylim([-0.5 0.5])
%     title(num2str(ii));
%     hold off
%     subplot(2,2,2);
%     plot(histgrad);
%     ylim([-0.5 0.5])
%     hold off
    
    nhistgrad =(histgrad./max(histgrad));
    zerointersectionlow = find(abs(nhistgrad)<0.1);
    zerointersectionhigh = find(nhistgrad>0.1,1,'first');
    zerointersectionlow(zerointersectionlow>zerointersectionhigh)=[];
    if isempty(zerointersectionlow)
        intersectz=min(zerointersectionhigh)-1;
    else
    intersectz = round(mean(zerointersectionlow));
    end
        
    
%     if ~(length(locs)<1) && ~(length(fpvec)<100)
    if ~(length(locs)<1) && ~(length(fpvec)<25)
        threshm = zeros(1,3);
        threshm(1) = Eedges(locs(1));
        threshm(2) = Eedges(intersectz);
        threshm(3) = mean([Eedges(locz(1)) Eedges(locs(1))]);
%         thresh = min(threshm); %original
%         thresh = mean(threshm);
        thresh = max(threshm);
        threshmat(ii) = thresh;
    elseif ~(length(locs)<1) && ~(length(fpvec)<10)
        thresh = Eedges(locs(1)+2);
        threshmat(ii) = thresh;
    elseif ~(length(locs)<1) && ~(length(fpvec)<3)
        thresh = Eedges(locs(1)+4);
        threshmat(ii) = thresh;
    else
        threshmat(ii) = NaN;
    end
% 
% 
%     if sum(ii == 15:5:25)
%     figure(99)
%     subplot(2,2,1);plot(Eedges(1:end-1),Ncc); hold all
%     if ~isnan(threshmat(ii))
%     stem(thresh,ones(size(thresh)).*0.3);
%     end
%     end
end

    
a = find(isnan(threshmat)==1);
    if ~isempty(a)
              threshmat(a) = ones(size(a)).*nanmax(threshmat);
    end
threshmatz = smooth(threshmat);

dottedz = false(size(fishdotstack));    
    for ii = 1:size(fishdotstack,3)
        ff = fishdotstack(:,:,ii);
        dot = zeros(size(ff));
        dot(ff>threshmatz(ii)) = 1;
        dottedz(:,:,ii) = dot;
    end    
        

savemstack(dottedz,details,SAV)
end



function [details,fishdotstack] = makeimagestack(EXP,pvalue,chan,ksize)
cd (EXP)
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


  %%%just load up 2 images instead of all the images
filename = char(imagefilenames{1});

[aa,bb] = regexp(filename,'p[0-9]+');
details.pvalue = filename(aa:bb);

[cc,dd] = regexp(filename,'[0-9]hr');
details.tpoint = filename(cc:dd);

[ee,ff] = regexp(filename,'(off|low|medlow|med|high|0dot00|0dot01|0dot02|0dot03|0dot04|0dot07|2dot40)');
dose = filename(ee:ff);
    if strcmp(details.tpoint,'0hr')
    dose = '0dot00';
    end
details.dose = dose;

[gg,hh] = regexp(filename,'(594_snail|594_smad7|594_pai1|647_smad7|594_pmepa1|594_tieg|594_bhlhe40|647_pai1|647_wnt9a|647_ctgf|594_wnt9a|594_ctgf|647_snail)');
species = filename(gg:hh);
details.chanspecie = species;
details.species = species;
details.channelnumber = chan;

 %collect date information  
[~,jj] = regexp(filename,'flat_2');
details.date = filename(jj:jj+9);

disp(strcat(details.species,details.pvalue));

dims = details.dims;
fishdotstack = LaplacianOfGaussianStack(stack,dims,ksize); %function in code
stphre=1;


end


function  LoGstack = LaplacianOfGaussianStack(imgstack,dims,ksize)
        LoGstack = zeros(dims(1),dims(2),dims(3));
        for i = 1:size(imgstack,3)
        LoGstack(:,:,i) = logMasked(imgstack(:,:,i),ksize);
        end
        
end


function savemstack(mstack,details,SAVdir)
cd(SAVdir)
save(strcat('mstack_',details.species,'_',details.dose,'_',details.tpoint,'_',details.channelnumber,'_',details.pvalue),'mstack');
cd ..
end

function savemrna(CCmrna,details,thresh,SAVdir)
cd(SAVdir)

save(strcat('details_',details.species,'_',details.dose,'_',details.tpoint,'_',details.channelnumber,'_',details.pvalue),'details');
save(strcat('CCmrna_',details.species,'_',details.dose,'_',details.tpoint,'_',details.channelnumber,'_',details.pvalue),'CCmrna');
save(strcat('thresh_',details.species,'_',details.dose,'_',details.tpoint,'_',details.channelnumber,'_',details.pvalue),'thresh');
cd ..
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