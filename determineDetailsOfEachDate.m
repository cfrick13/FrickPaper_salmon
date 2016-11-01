function determineDetailsOfEachDate
mfile = mfilename('fullpath');
[~,b] = regexp(mfile,'FrickPaperData');
mfiledir = mfile(1:b+1);
parentdir = mfiledir;

% parentdir = 'F:\2016_04_18 backup of new comp\pictures\Frick';

subdirname = '*FISHareaNewThresh5*';
uniqueconditions = struct();

cd(parentdir)
dirlist = dir('*201*smFISH');
dirnames = {dirlist.name};
%     for dnum = 1:length(dirnames)
    for dnum = 1:length(dirnames)
        cd(parentdir)
        dname = dirnames{dnum};
        cd(dname)
        
        subdirlist = dir(subdirname);
        if ~isempty(subdirlist)
            subdirnames = {subdirlist.name};
            subdname = subdirnames{1};
            cd(subdname)
            filelist = dir('*CCcells*');
            filename = {filelist.name};
                [~,~,~,d] = regexp(filename,'(594_snail|594_smad7|594_pai1|647_snail|647_pai1|647_smad7|647_wnt9a|647_ctgf|594_wnt9a|594_ctgf|647_smad7|647_pai1|594_pmepa1|594_tieg|594_bhlhe40)');
                probe = cellfun(@(x) x{1},d,'UniformOutput',0);
                [~,~,~,d] = regexp(filename,'(off|low|medlow|med|high|0dot00|0dot02|0dot03|0dot04|0dot07|2dot40)');
                dose = cellfun(@(x) x{1},d,'UniformOutput',0);
                [~,~,~,d] = regexp(filename,'(off|low|medlow|med|high|[0-9]dot[0-9]+)');
                dosealt = cellfun(@(x) x{1},d,'UniformOutput',0);
                [~,~,~,d] = regexp(filename,'[0-9]+hr');
                tpoints = cellfun(@(x) x{1},d,'UniformOutput',0);
                [~,~,~,d] = regexp(filename,'p[0-9]+');
                positions = cellfun(@(x) x{1},d,'UniformOutput',0);
                

        
            fulldetails = cell(1,length(positions));
            positiondetails = cell(1,length(positions));
            for i=1:length(positions)
               fulldetails{i} = strcat(probe{i},'-',dosealt{i},'-',tpoints{i},'-',positions{i});
               positiondetails{i} = strcat(probe{i},'-',dosealt{i},'-',tpoints{i});
            end

            details.probes = unique(probe);
            details.doses = unique(dose);
            details.tpoints = unique(tpoints);
            details.positions = unique(positions);
            details.fulldetails = fulldetails;
            uniqueconds = unique(positiondetails);
            for j = 1:length(uniqueconds)
                fieldname = uniqueconds{j};
                [~,~,~,d] = regexp(fulldetails,fieldname);
                pidx = cellfun(@isempty,d,'UniformOutput',1);
                pnames = positions(~pidx);
                fnames = fieldnames(uniqueconditions);
                    if isempty(fnames)
                        a=0;
                    else
                        a = length(uniqueconditions);
                    end
                uniqueconditions(a+1).date = dname;
                uniqueconditions(a+1).conditions = fieldname;
                uniqueconditions(a+1).pnames = horzcat(pnames{:});
            end
            
        else
            probe=[];
            dose=[];
            dosealt=[];
            tpoints=[];
            positions=[];
            
        end
        stophere=1;
    end

    
    snailspecific = determineSpecific(uniqueconditions,'snail');
    ctgfspecific = determineSpecific(uniqueconditions,'ctgf');
        stophere=1;


end


function structOut = determineSpecific(uniqueconditions,desireddetail)
conditions = {uniqueconditions(:).conditions}';
[~,~,~,d] = regexp(conditions,desireddetail);
didx = cellfun(@isempty,d,'UniformOutput',1);
structOut = uniqueconditions(~didx);
end
