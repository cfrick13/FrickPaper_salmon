function smfishmoveflatfield(DateB)
for B = DateB
    Date = B{1};
mfile = mfilename('fullpath');
[~,b] = regexp(mfile,'FrickPaperData');
mfiledir = mfile(1:b+1);
parentdir = mfiledir;
A = strcat(parentdir,Date,' smFISH\');
cd(A)
% [~,~,raw] = xlsread('smFISHdetailz.xlsx');

% for ii = 1:size(raw,1)-1 %if experiment number is numeric and not a string
%     r = raw{ii+1,1};
%     rnum = num2str(r);
%     raw{ii+1,1} = rnum;
% end

mkdir('FLATFIELDold')
cd('FLATFIELDold')
B = pwd;

cd ..
cd('FLATFIELD') %EXPORT ZEN FILES INTO THIS FOLDER
% fishfile = dir('*experiment*32*tif*');
fishfile = dir('*tif*');
filelistall = {fishfile.name};

datz = {fishfile.date};
[~,~,~,d] = regexp(datz,'\w+-[0-9][0-9][0-9][0-9]');
dmat = cellfun(@(x) x{1},d,'UniformOutput',0);
udmat = unique(dmat);
[a,b,c,d] = regexp(dmat,'2016');
ddmat = cellfun(@isempty,d,'UniformOutput',1);
filelistold = filelistall(ddmat);

parfor i=1:length(filelistold)

    fishfilename = char(filelistold{i});
    newfishfilename = fishfilename;
    disp(newfishfilename)
%     copyfile(fishfilename,strcat(B,'\',char(newfishfilename)));
    movefile(fishfilename,strcat(B,'\',char(newfishfilename)));

end
end
end