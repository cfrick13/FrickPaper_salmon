
TFCcat = struct();
TABScat = struct();
for Datesz = {'2014_03_15'}
        
        AllCellTimeTracesMakerMacWinX(Datesz)
        CellzFileReadExperimentStructureMacX
        cd('D:\Users\zeiss\Documents\MATLAB\AllImagingDataCompiled')
        T = load('tracesstructfc.mat', 'Ttracesstruct');
        Ttracesstruct = T.Ttracesstruct;
        T = load('tracesstructabs.mat', 'tracesstruct');
        tracesstruct = T.tracesstruct;
        TFCcat.(strcat('i',char(Datesz),'i')) = Ttracesstruct.i2dot40i;
        TABScat.(strcat('i',char(Datesz),'i')) = tracesstruct.i2dot40i;
        cd('D:\Users\zeiss\Documents\MATLAB\')
        
        
end