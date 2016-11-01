function smfishrenameexports(DateB)
for B = DateB
    Date = B{1};
    
mfile = mfilename('fullpath');
[~,b] = regexp(mfile,'FrickPaperData');
mfiledir = mfile(1:b+1);
parentdir = mfiledir;
A = parentdir;

cd(A)
[~,~,raw] = xlsread('smFISHdetailz.xlsx');

for ii = 1:size(raw,1)-1 %if experiment number is numeric and not a string
    r = raw{ii+1,1};
    rnum = num2str(r);
    raw{ii+1,1} = rnum;
end

mkdir('EXPORTS')
cd('EXPORTS')
B = pwd;

cd ..
cd('EXPORTSDIRECT') %EXPORT ZEN FILES INTO THIS FOLDER
% fishfile = dir('*experiment*32*tif*');
fishfile = dir('*tif*');
filelist = {fishfile.name};
% for i=1:length(filelist)
parfor i=1:length(filelist)
    fishfilename = char(filelist{i});
    [~,~,~,experiment] = regexp(fishfilename,'Experiment-[0-9][0-9]+');
    [~,~,~,experimentz] = regexp(char(experiment),'[0-9]+');
    experimentzz = char(experimentz);
%     experimentzz = str2num(experiment);
    [~,~,~,svalue] = regexp(fishfilename,'s[0-9]+');
    svalue = char(svalue);
    [~,~,~,channel] = regexp(fishfilename,'(TL DIC|Alexa Fluor 594|Alexa Fluor 647)');
    channel = char(channel);
    [~,~,~,zvalue] = regexp(fishfilename,'z[0-9][0-9]');
    zvalue = char(zvalue);
    


% disp(svalue)
% disp(experimentzz)

    expidx = strcmp(raw,experimentzz);
    [exprow,~] = find(expidx==1);
    svidx = strcmp(raw,svalue);
    [svrow,~] = find(svidx==1);
    expdetailsidx = ismember(svrow,exprow);
    expdetailsrow = svrow(expdetailsidx);
    disp(expdetailsrow)
    
    ExpDate = char(raw(expdetailsrow,3));
%     CellType = char(raw(expdetailsrow,4));
    FishProbe = char(raw(expdetailsrow,5));
    Row = char(raw(expdetailsrow,6));
    Well = char(raw(expdetailsrow,7));
    Dose = char(raw(expdetailsrow,8));
    TimePoint = char(raw(expdetailsrow,9));
    Pvalue = char(raw(expdetailsrow,11));
%     Fluorophore = char(raw(expdetailsrow,12));
    Reference = raw(expdetailsrow,10);
        
        if strcmp(channel,'Alexa Fluor 647')
            FishProbe = char(raw(expdetailsrow,13));
        end
        
%     disp(Reference)
        if isnan(Reference{1})
        newfishfilename = strcat(ExpDate,' smFISH ',{' '},FishProbe,{' '},Row,{' '},Well,{' '},'dose',{' '},Dose,{' '},TimePoint,{' '},Pvalue,'-',channel,'_',zvalue,'_ORG.tif');
        else
        newfishfilename = strcat(ExpDate,' smFISH ',{' '},FishProbe,{' '},Row,{' '},Well,{' '},'dose',{' '},Dose,{' '},TimePoint,{' '},'reference',{' '},Pvalue,'-',channel,'_',zvalue,'_ORG.tif');        
        end

%     copyfile(fishfilename,strcat(B,'\',char(newfishfilename)));
    disp(char(newfishfilename))
end
end
end