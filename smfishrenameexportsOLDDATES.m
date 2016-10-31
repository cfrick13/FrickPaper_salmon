function smfishrenameexportsOLDDATES(DateB)
for B = DateB
    Date = B{1};
A = strcat('D:\Users\zeiss\Pictures\Frick\',Date,' smFISH\');
cd(A)
% [~,~,raw] = xlsread('smFISHdetailz.xlsx');

% for ii = 1:size(raw,1)-1 %if experiment number is numeric and not a string
%     r = raw{ii+1,1};
%     rnum = num2str(r);
%     raw{ii+1,1} = rnum;
% end

mkdir('EXPORTS')
cd('EXPORTS')
B = pwd;

cd ..
cd('EXPORTSDIRECT') %EXPORT ZEN FILES INTO THIS FOLDER
% fishfile = dir('*experiment*32*tif*');
fishfile = dir('*tif*');
filelist = {fishfile.name};
parfor i=1:length(filelist)
% for i=1:length(filelist)
    fishfilename = char(filelist{i});
%     disp(fishfilename)
%     [~,~,~,experiment] = regexp(fishfilename,'Experiment-[0-9][0-9]+');
%     [~,~,~,experimentz] = regexp(char(experiment),'[0-9]+');
%     experimentzz = char(experimentz);
%     experimentzz = str2num(experiment);
    [~,~,~,svalue] = regexp(fishfilename,'p[0-9]+');
    svalue = char(svalue);
    svalue(1) = 'p';
    
    if isempty(svalue)
        [~,~,~,svalue] = regexp(fishfilename,'s[0-9]+');
        svalue = char(svalue);
    end
    
    [~,~,~,channel] = regexp(fishfilename,'(TL DIC|Alexa Fluor 594|Alexa Fluor 647|EGFP|EBFP2)');
    channel = char(channel);
    [~,~,~,zvalue] = regexp(fishfilename,'z[0-9][0-9]');
    zvalue = char(zvalue);
    
%     disp(channel)
    channelnum = channel(end-2:end);

    
    [~,~,~,d] = regexp(fishfilename,'(594_snail|594_smad7|594_pai1|647_snail|647_wnt9a|647_ctgf|594_wnt9a|594_ctgf|647_smad7|647_pai1|594_bhlhe40)');
    dcell = cellfun(@(x) x,d,'UniformOutput',0);

    if ~strcmp(channelnum,'594') & ~strcmp(channelnum,'647'); %this is if DIC is the channel
        channelnum = '594';
    end
    
    [~,~,~,d] = regexp(dcell,channelnum);
    didx = cellfun(@isempty,d,'UniformOutput',1);
    FishProbe = char(dcell(~didx));
    
    %determine row and well for 2015_01_19
    if strcmp(Date,'2015_01_19')
        snum = num2str(svalue(end-1:end));
        if snum>27
            snum=snum-27;
        end
        if snum <3
            Row = 'row2';
            Well = 'well1';
        elseif snum<5
            Row = 'row2';
            Well = 'well2';
        elseif snum<7
            Row = 'row2';
            Well = 'well3';
        elseif snum<9
            Row = 'row2';
            Well = 'well4';  
        elseif snum<11
            Row = 'row2';
            Well = 'well5';  
        elseif snum<13
            Row = 'row2';
            Well = 'well6';  
        elseif snum<15
            Row = 'row3';
            Well = 'well6';
        elseif snum<17
            Row = 'row3';
            Well = 'well5';
        elseif snum<19
            Row = 'row3';
            Well = 'well4';
        elseif snum<21
            Row = 'row3';
            Well = 'well3';
        elseif snum<23
            Row = 'row3';
            Well = 'well2';
        elseif snum<28
            Row = 'row3';
            Well = 'well1';
        elseif snum<19
            Row = 'row3';
            Well = 'well1';
        else
            Row = 'row4';
            Well = 'well7';
        end
        
    elseif strcmp(Date,'2015_01_29')
       if ~isempty(regexp(fishfilename,'smad7'))
            Row = 'row1';
            [~,~,~,d] = regexp(fishfilename,'well[0-9]+');
            Well = d{1};
       else
            [~,~,~,d] = regexp(fishfilename,'row[0-9]+');
            Row = d{1};
            [~,~,~,d] = regexp(fishfilename,'well[0-9]+');
            Well = d{1};
       end
    else
        [~,~,~,d] = regexp(fishfilename,'row[0-9]+');
        Row = d{1};
        [~,~,~,d] = regexp(fishfilename,'well[0-9]+');
        Well = d{1};
    end

    
    %dose for 2015_01_31
    if strcmp(Date,'2015_01_31')
        [~,~,~,d] = regexp(fishfilename,'(off|low|med|high|reference)');
        Dose = d{1};
        if strcmp(Dose,'off'); Dose = '0dot00';
        elseif strcmp(Dose,'low');Dose = '0dot04';
        elseif strcmp(Dose,'med');Dose = '0dot07';
        elseif strcmp(Dose,'high');Dose = '2dot40';
        elseif strcmp(Dose,'reference'); 
            Dose = '0dot00';
        end
        
        [~,~,~,d] = regexp(fishfilename,'([0-9]hr|reference)');
        if strcmp(d,'reference')
            TimePoint = '0hr';
        else
            TimePoint = d{1};
        end
        
        [~,~,~,d] = regexp(fishfilename,'(snail|smad7|pai1|reference)');
        dcell = cellfun(@(x) x,d,'UniformOutput',0);
        if strcmp(dcell,'reference')
            FishProbe = strcat(channelnum,'_snail');
        elseif strcmp(channelnum,'647')
            FishProbe = '647_pai1';
        else
            [~,~,~,d] = regexp(dcell,channelnum);
            species = char(dcell);
            FishProbe = strcat(channelnum,'_',species);
        end
        
    elseif strcmp(Date,'2015_01_29')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%reference has nothing
%0hr has no dose
%1hr has "med" or "high"
%all 647 is pai1
%some are "smad7" or "snail"

        [~,~,~,d] = regexp(fishfilename,'(off|low|med|high|reference)');
        if isempty(d)
            Dose='0dot00';
        else
            Dose = d{1};
        end
        if strcmp(Dose,'off'); Dose = '0dot00';
        elseif strcmp(Dose,'low');Dose = '0dot04';
        elseif strcmp(Dose,'med');Dose = '0dot07';
        elseif strcmp(Dose,'high');Dose = '2dot40';
        elseif strcmp(Dose,'reference'); 
            Dose = '0dot00';
        end
        
        [~,~,~,d] = regexp(fishfilename,'([0-9]hr|reference)');
        if strcmp(d,'reference')
            TimePoint = '0hr';
        elseif strcmp(d,'0hr')
            Dose = '0dot00';
            TimePoint = d{1};
        else
            TimePoint = d{1};
        end
        
        [~,~,~,d] = regexp(fishfilename,'(snail|smad7|pai1|reference)');
        dcell = cellfun(@(x) x,d,'UniformOutput',0);
        if strcmp(dcell,'reference')
            FishProbe = strcat(channelnum,'_snail');
        elseif strcmp(channelnum,'647')
            FishProbe = '647_pai1';
        else

            species = char(dcell);
            FishProbe = strcat(channelnum,'_',species);
        end
        
        
    elseif strcmp(Date,'2015_01_19')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%reference has nothing 
%0hr is "off"
%[0-9]hr has "med" or "high"
%all 647 is pai1
%some are "smad7" or "snail"
        [~,~,~,d] = regexp(fishfilename,'(off|low|med|high|reference)');
        Dose = d{1};
        if strcmp(Dose,'off'); Dose = '0dot00';
        elseif strcmp(Dose,'low');Dose = '0dot04';
        elseif strcmp(Dose,'med');Dose = '0dot07';
        elseif strcmp(Dose,'high');Dose = '2dot40';
        elseif strcmp(Dose,'reference'); 
            Dose = '0dot00';
        end
        
        [~,~,~,d] = regexp(fishfilename,'([0-9]hr|reference)');
        if strcmp(d,'reference')
            TimePoint = '0hr';
        else
            TimePoint = d{1};
        end
        
        [~,~,~,d] = regexp(fishfilename,'(snail|smad7|pai1|reference)');
        dcell = cellfun(@(x) x,d,'UniformOutput',0);
        if strcmp(dcell,'reference')
            FishProbe = strcat(channelnum,'_snail');
        elseif strcmp(channelnum,'647')
            FishProbe = '647_pai1';
        else
            [~,~,~,d] = regexp(dcell,channelnum);
            species = char(dcell);
            FishProbe = strcat(channelnum,'_',species);
        end    
        
        
        
    %dose for 2015_03_06
    elseif strcmp(Date,'2015_03_06')
        [~,~,~,d] = regexp(fishfilename,'(0dot00|0dot01|0dot02|0dot03|0dot04|0dot07|2dot40|off|reference)');
        Dose = d{1};
        if strcmp(Dose,'off'); Dose = '0dot00';
        elseif strcmp(Dose,'reference'); 
            Dose = '0dot00';
        end
        
        [~,~,~,d] = regexp(fishfilename,'([0-9]hr|reference)');
        if strcmp(d,'reference')
            TimePoint = '0hr';
        else
            TimePoint = d{1};
        end
        
        [~,~,~,d] = regexp(fishfilename,'(reference|594_snail|594_smad7|594_pai1|647_snail|647_wnt9a|647_ctgf|594_wnt9a|594_ctgf|647_smad7|647_pai1|594_bhlhe40)');
        dcell = cellfun(@(x) x,d,'UniformOutput',0);
        if strcmp(dcell,'reference')
            FishProbe = strcat(channelnum,'_snail');
        else
            [~,~,~,d] = regexp(dcell,channelnum);
            didx = cellfun(@isempty,d,'UniformOutput',1);
            FishProbe = char(dcell(~didx));
            stophere=1;
        end
    else
        [~,~,~,d] = regexp(fishfilename,'(0dot00|0dot01|0dot02|0dot03|0dot04|0dot07|2dot40)');
        Dose = d{1};
        
        [~,~,~,d] = regexp(fishfilename,'[0-9]hr');
        TimePoint = d{1};
    end
    

        
    Pvalue = svalue;
    ExpDate = Date;
    
    [~,~,~,d] = regexp(fishfilename,'reference');
    if ~isempty(d)
        Reference = d(1);
    else
        Reference{1} = NaN;
    end

        if isnan(Reference{1})
        newfishfilename = strcat(ExpDate,' smFISH ',{' '},FishProbe,{' '},Row,{' '},Well,{' '},'dose',{' '},Dose,{' '},TimePoint,{' '},Pvalue,'-',channel,'_',zvalue,'_ORG.tif');
%             if strcmp(zvalue,'z10')
%             disp(char(newfishfilename))
%             end
            disp(char(newfishfilename))
        else
        newfishfilename = strcat(ExpDate,' smFISH ',{' '},FishProbe,{' '},Row,{' '},Well,{' '},'dose',{' '},Dose,{' '},TimePoint,{' '},'reference',{' '},Pvalue,'-',channel,'_',zvalue,'_ORG.tif');        
%             if strcmp(zvalue,'z10')
%             disp(char(newfishfilename))
%             end
            disp(char(newfishfilename))
        end
%     disp(char(newfishfilename))
    copyfile(fishfilename,strcat(B,'\',char(newfishfilename)));

end
end
end