director = '\FISHareaNewTHRESH';
% director = '\FISHareaNewMEDIAN';

% Date = {'2016_02_09','2016_02_11','2016_02_19'};
% Date = {'2016_02_09','2016_02_11'};

% Date = {'2015_08_31','2015_09_03','2015_12_11','2015_12_15','2015_12_19','2016_01_25','2016_02_09','2016_02_11'};
Date = {'2015_01_19','2015_01_29','2015_01_31','2015_03_06','2015_03_25','2015_03_31','2015_04_01','2015_08_31','2015_09_03','2015_12_11','2015_12_15','2015_12_19','2016_01_25','2016_02_09','2016_02_11'};
% Date = {'2015_03_06'};
% Date = {'2015_01_31','2015_03_06','2015_03_25','2015_03_31','2015_04_01'};
% Date = {'2015_01_31','2015_03_06'};
% Date = {'2015_12_11','2015_12_15','2015_12_19','2016_01_25','2016_02_09','2016_02_11'};
% Date = {'2015_12_11','2015_12_15','2015_12_19','2016_01_25'};

% Date = {'2015_01_19','2015_01_29'};
smfishrenameexportsOLDDATES(Date) % renameFISHimages
smfishmoveflatfield(Date)
smfishrenameexports(Date) % renameFISHimages
smfishflatbkgPARALLELnew(Date)
% % % LetsGoFishPar
% % % LetsGoFishParATTEMPT(Date)
LetsGoFishParATTEMPTrefined(Date,director)

%     %check these images with uiMRNAimagesMRNAmstackWin
%     uiMRNAimagesMRNAdottedWin(director)
    mstackTOdottedWin(Date,director)
getSavedFocuses(Date)
% toughsegmentationforIMMUNOforFISH(Date)

% uicorrectSegmentedImagesWin
LetsGoFishPlot(Date,director)


% uiPositionsNew
writeEmmiesNew(director)
AllCellTimeTracesMakerMacWin
CellzFileReadExperimentStructureMac
plotCORRELATIONSfromStructureBOXPLOTcompare


focusEmmies
testfocusFISH
plotMRNAnew(Date,director)

%%




%% Do not edit
% smfishrenameexports(Date) % renameFISHimages
% smfishflatbkgPARALLEL(Date)
% % LetsGoFishPar
% LetsGoFishParATTEMPT(Date)
%     %check these images with uiMRNAimagesMRNAmstackWin
%     %uiMRNAimagesMRNAmstackWin
% getSavedFocuses(Date)
% toughsegmentationforIMMUNOforFISH(Date)
% uicorrectSegmentedImagesWin
% plotMRNAnew

% uiPositionsNew
% LetsGoFishPlot
% writeEmmiesNew



%% Do not edit old
% smfishrenameexports(Date) % renameFISHimages
% smfishflatbkgPARALLEL(Date)
% % LetsGoFishPar
% LetsGoFishParATTEMPT(Date)
%     %check these images with uiMRNAimagesMRNAmstackWin
%     %uiMRNAimagesMRNAmstackWin
% getSavedFocuses(Date)
% toughsegmentationforIMMUNOforFISH(Date)
% 
% uicorrectSegmentedImagesWin
% uiPositionsNew
% LetsGoFishPlot
% writeEmmiesNew
% 
% focusEmmies
% testfocusFISH