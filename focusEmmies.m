function focusEmmies

CountsString = 'mALLFOCUS.mat';
emmiesfiles = 'mALL.mat';
A = 'D:\Users\zeiss\Pictures\Frick\';
% for BB = {'2015_08_31 smFISH','2015_09_03 smFISH'};
for BB = {'2015_01_15 smFISH','2015_01_19 smFISH','2015_01_29 smFISH','2015_01_31 smFISH','2015_03_06 smFISH','2015_03_25 smFISH','2015_03_31 smFISH','2015_04_01 smFISH','2015_08_31 smFISH','2015_09_03 smFISH','2015_08_31 smFISH','2015_09_03 smFISH','2015_12_15 smFISH','2015_12_19 smFISH','2016_01_25 smFISH','2016_02_09 smFISH','2016_02_19 smFISH'};
    B=BB{1};

    
    
    
    
cd('D:\Users\zeiss\Documents\MATLAB\AllImagingDataCompiled')

mal = load(emmiesfiles);
fnamesmal = fieldnames(mal);
dateofexp = B(1:10);    
    
cd (strcat(A,B))
filelist = dir('*infocuspositions.xlsx');
if ~isempty(filelist)
    filename = char(filelist.name);
    [num,txt,raw] = xlsread(filename);
    pvalues = raw;
else
    %set all positions as in focus
    for i=1:40; pstr = num2str(i); pchar = 'p00'; pchar(end-(length(pstr)-1):end)=pstr; pvalues{i}=pchar;
    end
end



for ubi=1:length(pvalues)
            if ubi==1
            pvaluelist = strcat('(',pvalues{ubi},'|');
            elseif ubi==length(pvalues)
                pvaluelist = horzcat(pvaluelist,strcat(pvalues{ubi},')'));
            else
                pvaluelist = horzcat(pvaluelist,strcat(pvalues{ubi},'|'));
            end
end

[~,~,~,d,~] = regexp(fnamesmal,pvaluelist);
didxpvalue = cellfun(@isempty,d,'UniformOutput',1);

[~,~,~,d,~] = regexp(fnamesmal,dateofexp);
didxdate = cellfun(@isempty,d,'UniformOutput',1);
mmm = ~didxpvalue & ~didxdate;
mnidx = find(mmm==1);
cd('D:\Users\zeiss\Documents\MATLAB\')

for nidx = mnidx'
    emmiesname = fnamesmal{nidx};
        if isempty(dir(CountsString)) %write original
            save(CountsString, '-struct','mal',emmiesname);
        else    %append
            save(CountsString, '-struct','mal',emmiesname,'-append');
        end
end
        
end
end
