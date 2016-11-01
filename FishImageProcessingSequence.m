director = '\FISHareaNewTHRESH';
% director = '\FISHareaNewMEDIAN';

% Date = {'2016_02_09','2016_02_11','2016_02_19'};
% Date = {'2016_02_09','2016_02_11'};
% Date = {'2015_03_06'};
% Date = {'2015_01_31','2015_03_06','2015_03_25','2015_03_31','2015_04_01'};
% Date = {'2015_01_31','2015_03_06'};
% Date = {'2015_12_11','2015_12_15','2015_12_19','2016_01_25','2016_02_09','2016_02_11'};
% Date = {'2015_12_11','2015_12_15','2015_12_19','2016_01_25'};
% Date = {'2015_08_31','2015_09_03','2015_12_11','2015_12_15','2015_12_19','2016_01_25','2016_02_09','2016_02_11'};

% %HELPFUL FILES THAT COULD BE NECESSARY
% smfishrenameexportsOLDDATES(Date) % renameFISHimages
% smfishmoveflatfield(Date)

% %FIRST LOOP
% Date = {'2015_01_19','2015_01_29','2015_01_31','2015_03_06','2015_03_25','2015_03_31','2015_04_01','2015_08_31','2015_09_03','2015_12_11','2015_12_15','2015_12_19','2016_01_25','2016_02_09','2016_02_11'};
% smfishrenameexports(Date) % renameFISHimages
% smfishflatbkgPARALLELnew(Date)
% LetsGoFishParATTEMPTrefined(Date,director)
% 
% 
% %CHECK FILES TO SEE SUCCESS OF FIRST LOOP
% % % % %check these images with 
% % mstackTOdottedWin(Date,director)
% % uiMRNAimagesMRNAdottedWin(director)
% getSavedFocuses(Date) %might not be necessary actually
% 
% %PERFORM SEGMENTATION based on FISH fluorescence
% toughsegmentationforIMMUNOforFISH(Date)
% uicorrectSegmentedImagesWin

%SECOND LOOP
% Date = {'2015_01_19','2015_01_29','2015_01_31','2015_03_06','2015_03_25','2015_04_01','2015_12_11','2016_01_25'};%snail
% Date = {'2015_04_01'};%snail
Date = {'2016_01_25'};
% Date = {'2015_01_19','2015_01_29','2015_01_31','2015_03_06','2015_03_25','2015_03_31','2015_04_01','2015_08_31','2015_09_03','2015_12_11','2015_12_15','2015_12_19','2016_01_25','2016_02_09','2016_02_11'};
% LetsGoFishPlot(Date,director)


% %ASSIGN CELL NUMBERS BASED ON LIVE CELL IMAGING EXPERIMENTS
% uiPositionsNew
% 
% %COUNT MRNA IN CELLS
writeEmmiesNew(director)
% 
% % %ASSIGN MRNA COUNTS TO LIVE CELL IMAGING DATA and PLOT
% AllCellTimeTracesMakerMacWin
% CellzFileReadExperimentStructureMac
% plotCORRELATIONSfromStructureBOXPLOTcompare
% 
% %EXTRA FILES
focusEmmies
% testfocusFISH
% plotMRNAnew(Date,director)
