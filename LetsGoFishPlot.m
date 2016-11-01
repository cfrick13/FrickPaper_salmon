function LetsGoFishPlot(Date,director)

mfile = mfilename('fullpath');
[~,b] = regexp(mfile,'FrickPaperData');
mfiledir = mfile(1:b+1);
parentdir = mfiledir;
A = parentdir;
% for BB = {'2015_01_15 smFISH','2015_01_19 smFISH','2015_01_29 smFISH','2015_01_31 smFISH','2015_03_06 smFISH','2015_03_25 smFISH','2015_03_31 smFISH','2015_04_01 smFISH'};
%     for BB = {'2015_01_29 smFISH','2015_01_30 smFISH','2015_01_31 smFISH','2015_03_06 smFISH'};
% for BB = {'2015_01_15 smFISH','2015_01_19 smFISH','2015_01_29 smFISH','2015_01_31 smFISH','2015_03_06 smFISH','2015_03_25 smFISH','2015_03_31 smFISH','2015_04_01 smFISH','2015_05_14 smFISH'}
% for BB = {'2015_01_15 smFISH','2015_01_19 smFISH','2015_01_29 smFISH','2015_01_30 smFISH','2015_01_31 smFISH','2015_03_06 smFISH'}
% for BB = {'2015_01_15 smFISH','2015_01_19 smFISH','2015_01_29 smFISH','2015_01_31 smFISH','2015_03_06 smFISH','2015_03_25 smFISH','2015_03_31 smFISH','2015_04_01 smFISH','2015_05_14 smFISH','2015_07_02 smFISH','2015_07_10 smFISH'};
%      for BB = {'2015_12_15 smFISH','2015_12_19 smFISH'};
%     for BB = {'2015_08_31 smFISH'};
%        for BB = {'2015_01_30 smFISH'}
    for BB = Date
    
    


    B = BB{1};
    B = strcat(B,' smFISH');
    D = '\FLATFIELD';EE = '\autoseg';
    F = director;

        for ksize = 5;
        EXP = strcat(A,B,F,num2str(ksize));
        SAV = strcat(A,B,F,num2str(ksize));
        SEG = strcat(A,B,EE);
        FLAT = strcat(A,B,D);

        cd(EXP)
        filelist = dir(strcat('*thresh*.mat'));
        CHANS = findNumberOfVarsInList(filelist,'(594|647)');
%         CHANS = findNumberOfVarsInList(filelist,'(647)');
            for channels = CHANS
%             for channels = {'594','647'};                
            channel = char(channels);
            cd(SEG)
            filelist = dir(strcat('*sm*.mat'));
                if isempty(filelist)
                    filelist = dir(strcat('*C2*.mat'));
                end
            PVALUES = findNumberOfVarsInList(filelist,'p[0-9]+');
            
           
            disp(B)
                parfor cycle = 1:length(PVALUES)
%                 for cycle = 1:length(PVALUES)
                pvalue = char(PVALUES(cycle));
                cd(EXP)
                MRNAfilelist = dir(strcat('CCmrna*',channel,'*',pvalue,'.mat'));
                filename = char({MRNAfilelist.name});
                details = setDetails(filename,channel,B);
                CCmrna = loadupfile('CCmrna',details); %loads CCmrna
                %%%load up the segLogical
                cd(SEG)
                SEGfilelist = dir(strcat('*sm*',pvalue,'_*594.mat'));
                    if isempty(SEGfilelist)
                        SEGfilelist = dir(strcat('*sm_',pvalue,'.mat'));
                    end
                SEGfilename = char({SEGfilelist.name});
                ImageSize = CCmrna.ImageSize;
                    if ~isempty(dir(SEGfilename))
                        CCcells = loadUpSegLogical(ImageSize,SEGfilename,'segLogical');
                        cd(EXP)
                        CCfoplot = countMrna(CCcells,CCmrna,details,SAV,ksize,FLAT);
                    else
                        disp(strcat('no seg file',pvalue))
                    end
                end
            end
        end
    end
end



function CCcells = loadUpSegLogical(ImageSize,SEGfilename,str)
load(SEGfilename,str);
segmentedcellstack = zeros(ImageSize);
    for i = 1:ImageSize(3)
        segmentedcellstack(:,:,i) = segLogical;
    end
CCcells = bwconncomp((segmentedcellstack));
end

function details = setDetails(filename,channel,B)
[aa,bb] = regexp(filename,'p[0-9]+');
details.pvalue = filename(aa:bb);
[cc,dd] = regexp(filename,'[0-9]hr');
details.tpoint = filename(cc:dd);
[ee,ff] = regexp(filename,'(off|low|medlow|med|high|0dot00|0dot02|0dot03|0dot04|0dot07|2dot40)');
dosez = filename(ee:ff);
[gg,hh] = regexp(filename,'(smad7|snail|row1|row2|594_snail|594_smad7|594_pai1|647_snail|647_pai1|647_smad7|647_wnt9a|647_ctgf|594_wnt9a|594_ctgf|647_smad7|647_pai1|594_pmepa1|594_tieg|594_bhlhe40)');
details.species = filename(gg:hh);
details.channelnumber = channel;
details.date = B(1:11);
    if strcmp(details.tpoint,'0hr')
        dose = '0dot00';
    else
        dose=dosez;
    end
details.dose = dose;
details.rownumb = details.species;
end

function CCcells = countMrna(CCcells,CCmrna,details,SAVdir,~,~)
%     if strcmp(details.channelnumber,'594')
%         if ksize == 5;
% %         arealimit = 3;
%         arealimit = 1;
%         elseif ksize == 7;
%         arealimit = 3;
%         end
%     elseif strcmp(details.channelnumber,'647')
%         if ksize == 5;
%         arealimit = 1;
%         elseif ksize == 7;
%         arealimit = 3;
%         end 
%     end
arealimit =1;

%%%%%%%%%%

%%%%%%%%%

mRNAinCell = zeros(CCcells.NumObjects,1);
areaOfCell = zeros(CCcells.NumObjects,1);
centroid = cell(CCcells.NumObjects,1);
% statsMRNA = regionprops(CCmrna,'Area');
% statsCELLS = regionprops(CCcells,'centroid');
statsCELLS = struct();
PXX =  CCcells.PixelIdxList;
    for j = 1:length(PXX)
        px = PXX{j};
        [row,col,zs] = ind2sub(CCcells.ImageSize,px(1)); %x and y come out reverse of S.Centroid
        statsCELLS(j).Centroid = [col,row,zs];
    end
    

CellsWithColor = zeros(CCcells.ImageSize);
    for qq = 1:CCcells.NumObjects
       cellshape = CCcells.PixelIdxList{qq};
       areaOfCell(qq) =  numel(cellshape);
       cent = statsCELLS(qq).Centroid;
       centroid{qq} = [cent(1) cent(2)];
       CellsWithColor(cellshape) = qq;
    end

    

    px = CCmrna.PixelIdxList;
    pxl = arrayfun(@(x) length(x{1}),px,'UniformOutput',1);
    pxz = pxl>1;
    pxlist = px(pxz);
        CCmrna.NumObjects = length(pxlist);
        CCmrna.PixelIdxList = pxlist;


arealimit=1;

%%%cycle through each mRNA dot to assign it to a cell. 
    if ~isempty(CCmrna.PixelIdxList)
        for q = 1:length(CCmrna.PixelIdxList)
            if length(CCmrna.PixelIdxList{q}) < arealimit
            else
            dot = CCmrna.PixelIdxList{q};
            [~,~,z] = ind2sub(CCcells.ImageSize,dot);
            zul = length(unique(z));
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %             %%%% Make the DOTTED stack that circles mrna for count validation
% % %             [x,y,z] = ind2sub(CCcells.ImageSize,dot);
% % %             frankx{q} = x;
% % %             franky{q} = y;
% % %             frankz{q} = z;
% % %             meanx = round(nanmean(x));
% % %             meany = round(nanmean(y));
% % %             meanz = round(nanmean(z));
% % %                 if meanx==2048 || meanx ==1 || meany==1 || meany==2048
% % %                     dotted(meanx:meanx,meany:meany,meanz)=1;
% % % %                     disp('edge')
% % %                 else
% % %                     dotted(meanx:meanx,meany:meany,meanz)=1;
% % %                 end
% % %                 
% % %             %%%% Make the DOTTED stack that circles mrna for count validation
% % %             [x,y,z] = ind2sub(CCcells.ImageSize,dot(1));
% % %             meanx = x(1);
% % %             meany = y(1);
% % %             meanz = z(1);
% % %                 if meanx==2048 || meanx ==1 || meany==1 || meany==2048
% % %                     dotted(meanx:meanx,meany:meany,meanz)=1;
% % % %                     disp('edge')
% % %                 else
% % %                     dotted(meanx:meanx,meany:meany,meanz)=1;
% % %                 end    
% % %                 
% % %                 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            qq = round(mean(mean(mean(CellsWithColor(dot)))));
                if qq==0
                else
                    if zul<6
                        mRNAinCell(qq) = mRNAinCell(qq)+1;
                    elseif zul>19
                        mRNAinCell(qq) = mRNAinCell(qq)+0;
                    else
                        mRNAinCell(qq) = mRNAinCell(qq)+round(zul./3);
                    end
                    
                end
            end
        end
    end
    
    


  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %     % % make dotted images have circles around counted mrna
% % % %         if strcmp(details.species,'594_bhlhe40')
% % % % %         elseif strcmp(details.dose,'high')
% % % %         SE = strel('disk',4);
% % % %         doots = imdilate(dotted,SE);
% % % %         dooter = bwperim(doots,8);
% % % % 
% % % %         %write DOTTED images into stack
% % % %         cd ..
% % % %             for i=1:size(dotted,3)
% % % %                 if i==1
% % % %             imwrite(uint16(dooter(:,:,i)),strcat('DOTTED_',details.species,'_',details.dose,'_',details.tpoint,'_',details.channelnumber,'_',details.pvalue,'.TIFF'));
% % % %                 else
% % % %             imwrite(uint16(dooter(:,:,i)),strcat('DOTTED_',details.species,'_',details.dose,'_',details.tpoint,'_',details.channelnumber,'_',details.pvalue,'.TIFF'),'writemode','append');
% % % %             %imwrite(rgb2gray(a), 'myFile.TIFF', 'writemode', 'append')
% % % %             disp('dotted write')
% % % %                 end
% % % %             end
% % % % 
% % % %         cd(FLAT)
% % % %         fs = strcat('*flat*',details.pvalue,'-*',details.channelnumber,'*');
% % % %         filelist = dir(fs);
% % % %         i=1;
% % % %             for cfile = {filelist.name}
% % % %                 cd(FLAT)
% % % %                 img = imread(char(cfile));
% % % %                 cd ..
% % % %                 if i==1
% % % %                     imwrite(uint16(img),strcat('IMG_',details.species,'_',details.dose,'_',details.tpoint,'_',details.channelnumber,'_',details.pvalue,'.TIFF'));
% % % %                     i=2;
% % % %                 else
% % % %                     imwrite(uint16(img),strcat('IMG_',details.species,'_',details.dose,'_',details.tpoint,'_',details.channelnumber,'_',details.pvalue,'.TIFF'),'writemode','append');
% % % %                 end
% % % %             end
% % % % 
% % % % 
% % % %         cd autoseg\
% % % %         seglogfilelist = dir(strcat('*flat*',details.pvalue,'_*.mat'));
% % % %         label = labelmatrix(CCcells);
% % % %         slname = char({seglogfilelist.name});
% % % %         A = load(slname);
% % % %         slimg = A.segLogical;%seglogical image
% % % % 
% % % % 
% % % %         figure(737)
% % % %         % subplot(1,2,1);
% % % %         % imagesc(slimg);
% % % %         % for i=1:length(A.cellLocations)
% % % %         % xy = A.cellLocations{i};
% % % %         % if ~isempty(xy)
% % % %         % x=xy(1);
% % % %         % y=xy(2);
% % % %         % text(x,y,A.cellIDs{i});hold on
% % % %         % % text(x+20,y+20,num2str(mRNAinCell(i)));hold on
% % % %         % end
% % % %         % end
% % % % 
% % % % 
% % % %         % subplot(1,2,2);
% % % %         imagesc(label(:,:,1));
% % % %             for i=1:length(centroid)
% % % %             xy = centroid{i};
% % % %                 if ~isempty(xy)
% % % %                 x=xy(1);
% % % %                 y=xy(2);
% % % %                 % text(x,y,A.cellIDs{i});hold on
% % % %                 text(x,y,num2str(mRNAinCell(i)));hold on
% % % %                 end
% % % %             end
% % % % 
% % % %         end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%determine concentration of mrna in cell and assign values to the structure
mRNAinCellconc = mRNAinCell./areaOfCell;
CCcells.(strcat('mRNA',details.channelnumber,'inCell')) = mRNAinCell;
CCcells.(strcat('mRNA',details.channelnumber,'inCellconc')) = mRNAinCellconc;
CCcells.('areaOfCell') = areaOfCell;
CCcells.('centroid') = centroid;

%save the structure and display it...
cd(SAVdir)
disp(strcat('CCcells_',details.species,'_',details.dose,'_',details.tpoint,'_',details.channelnumber,'_',details.pvalue));
save(strcat('CCcells_',details.species,'_',details.dose,'_',details.tpoint,'_',details.channelnumber,'_',details.pvalue),'CCcells');
end

function [cellLocation,segmentedIDs] = uniquecells(locals,CCcells)
CCcells=CCcells.CCcells;
segmentedIDs = zeros(size(locals));
mmm=1;
cellLocation = [];
    for e=1:length(CCcells.PixelIdxList)
    cell = CCcells.PixelIdxList{e};
        if length(cell)>5000
        [x,y,z] = ind2sub(size(locals),cell);
        segmentedIDs(cell) = mmm;
        meanx = round(mean(x));
        meany = round(mean(y));
        cellLocation{mmm} = [meanx,meany];
        mmm=mmm+1;
        end
    end
end

function details = countmrnaOLD(cellLocation,segmentedIDs,CCmrna,details,thresh,projfish,locals,SAVdir,B,ksize)

if strcmp(details.channelnumber,'594')
if ksize == 5;
arealimit = 3;
elseif ksize == 7;
    arealimit = 6;
end
elseif strcmp(details.channelnumber,'647')
   if ksize == 5;
arealimit = 1;
elseif ksize == 7;
    arealimit = 3;
   end 
end


CCmrna = CCmrna.CCmrna;
thresh = thresh.thresh;
mRNAinCell = zeros(length(cellLocation),1);
areaOfCell = zeros(length(cellLocation),1);
mRNAinCellconc = zeros(length(cellLocation),1);
lengths = zeros(length(cellLocation),1);
stats = regionprops(CCmrna,'Area');


dotted = zeros(size(locals));
if ~isempty(CCmrna.PixelIdxList)
for q = 1:length(CCmrna.PixelIdxList)
if stats(q).Area < arealimit
else
dot = CCmrna.PixelIdxList{q};
[x,y,z] = ind2sub(size(locals),dot);
meanx = round(nanmean(x));
meany = round(nanmean(y));
meanz = round(nanmean(z));
lengths(q) = length(z);
zfocus(q) = meanz;
dotlength(q) = length(dot);
a = mean(zfocus);
% dotted(x,y,z)=1;
dotted(meanx,meany,meanz)=1;
e = segmentedIDs(meanx,meany,meanz);
if e<1
else
mRNAinCell(e) = mRNAinCell(e)+1;
end
% 
%         figure(1),
%         if length(z) <10
%         hhh = plot3(x,y,z,'LineWidth',3,'Color',[1 1 0.5]./length(dot));hold all
%             clor = get(hhh,'Color');
%         plot3(meanx,meany,meanz,'LineStyle','none','MarkerEdgeColor','none','MarkerFaceColor',clor,'MarkerSize',5,'Marker','o');
%             zlim([0 25])
%             xlim([0 2024])
%             ylim([0 2024])
%         end
end
end
%         view ([90 -90]) %same orientation as final image
%         title(strcat('mRNA = ', num2str(length(CCmrna.PixelIdxList)),', thresh  ',num2str(thresh),'......',details.tpoint,'-',details.pvalue,'-',details.channelnumber));
%         cd(SAVdir)
%         filename = strcat('mRNAdots','_',details.rownumb,'_',details.tpoint,'_',details.channelnumber,'_',details.pvalue,'.fig');
%         savefig(filename)
%         close 1

for i=1:length(mRNAinCell)

OfCell = segmentedIDs(segmentedIDs==i);
areaOfCell(i) = numel(OfCell)./size(segmentedIDs,3);
end

mRNApercell = mean(mRNAinCell);
mRNAinCellconc = mRNAinCell./areaOfCell;

details.meanlength = mean(lengths);
a = mean(zfocus);
details.meanfocus = mean(zfocus);
details.mrnapercell = mRNApercell;
details.thresh = thresh;
details.dotlength = dotlength
% [numbers,bins] = hist(dotlength,[0:1:10]);hold all
% bar(bins,numbers);

%% plot mRNA only in segmented cells
% figure(2)
% cMap(1,1:3) = [0 0 0]; cMap(2:length(cellLocation)+2,1:3) = jet(length(cellLocation)+1);zslice=10;
% colorcells = ind2rgb((uint8(segmentedIDs(:,:,zslice))),cMap);
% cMap = gray(max(max(projfish)));
% 
% colordots = ind2rgb(projfish,cMap);
% h = imagesc(projfish);
% colormap('pink')
% hold on
% hh = image(colorcells);
% hold on
% set(hh,'AlphaData',0.3);
% hold off
%  set(gca,'YDir','normal');
% 
% 
% for e= 1:length(cellLocation)
% xy = cellLocation{e};
% ms = text(xy(2),xy(1),num2str(mRNAinCell(e)),'FontSize',18,'Color',[0.8 0.8 0.8]);
% mss = text(xy(2),xy(1)-50,num2str(mRNAinCellconc(e).*1000),'FontSize',8,'Color',[0.7 0.7 0.7]);
% title(strcat('mRNAperCell = ',num2str(mRNApercell),'......',details.dose,'_',details.tpoint,'-',details.pvalue,'-',details.channelnumber));
% end
% 
% 
% 
% 
% filename = strcat('CountsPerCell','_',details.rownumb,'_',details.dose,'_',details.tpoint,'_',details.channelnumber,'_',details.pvalue,'.fig');
% savefig(filename)
% close gcf

if strcmp(details.pvalue,'p10')
for i=1:size(dotted,3)
imwrite(dotted(:,:,i),strcat('smalldotted',num2str(i),'_',details.channelnumber,'_',details.pvalue,'.tif'),'Tiff');
end
end

cd(SAVdir)
filename = strcat(B,'_',details.rownumb,'_',details.dose,'_',details.tpoint,'_',details.channelnumber,'_',details.pvalue,'_',num2str(thresh),'.mat');
save(filename,'mRNAinCell','mRNAinCellconc','areaOfCell');




% 
%             figure(3)
%             [number,bincenters] = hist(lengths,[1:20]);
%             normnumber =sum(number);
%             fract = number./normnumber;
%             subplot(10,8,mm);bar(bincenters,fract);
%             ylim([0 0.8]);
%             mm=mm+1;
%             xlim([0 20]);
%             title(strcat('m/c=', num2str(mRNApercell),',t=',num2str(thresh),' , ',details.tpoint,'-',details.pvalue,'-',details.channelnumber));
%             filename = strcat('mRNAlengthHIST','_',details.tpoint,'-',details.pvalue,'-',details.channelnumber,'_thresh_',num2str(thresh),'.fig');
%             xlabel('mRNA length')
%             ylabel('number of mRNA');
%             savefig(filename)

end
end

function var = loadupfile(variable,details)
filelist = dir(strcat('*',variable,'*',details.tpoint,'*',details.channelnumber,'*',details.pvalue,'.mat'));
name = char({filelist.name});
A = load(name);
var = A.(variable);
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