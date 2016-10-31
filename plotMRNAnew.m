function plotMRNAnew(Date,director)
SpeciesToPlot = '594_snail';
COLORZ = {'k','r','y','g','b','m','c'};
SpeciesToCorrelate = {'ctgf','wnt9a'}; %for correlation!!!!

A = 'D:\Users\zeiss\Pictures\Frick\';
% for BB = {'2015_12_11','2015_12_15','2015_12_19','2016_01_25'};
% for BB = {'2015_12_15 smFISH'};
% for BB = {'2015_01_15 smFISH','2015_01_19 smFISH','2015_01_29 smFISH','2015_01_31 smFISH','2015_03_06 smFISH','2015_03_25 smFISH','2015_03_31 smFISH','2015_04_01 smFISH','2015_12_11 smFISH','2015_12_15 smFISH','2015_12_19 smFISH','2016_01_25 smFISH'};
% for BB = {'2016_02_09 smFISH','2016_02_11 smFISH'};

cycle=1;
for BB = Date
    
     B = BB{1};
    B = strcat(B,' smFISH');
    
    F = director;
        ksize = 5;
        EXP = strcat(A,B,F,num2str(ksize));
        cd(EXP)
    files = strcat('*CCcells*');
    filelist = dir(files);
    cfile = {filelist.name};

    %%%%%
%     M = runMRNA(cfile,B);
%     cd (strcat(A,B))
%     save('mrna.mat','M');
%     %
    cd (strcat(A,B))
    load('mrna.mat','M');
    %%%%%
    
    MM{cycle} = M;
    cycle=cycle+1;
end


stophere=1;
for i=1:length(MM)
    M = MM{i};
        if i==1
        mmrna = M;
        else
        mmrna = horzcat(M,mmrna);
        end
end

M=mmrna;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%plot individual scenes timecourse
for i=1:length(M)
    species = M(i).species;
    dose = M(i).dose;
        if strcmp(dose,'0dot00');d=1;
        elseif strcmp(dose,'0dot01');d=2;
        elseif strcmp(dose,'0dot02');d=3;    
        elseif strcmp(dose,'0dot03');d=4; 
        elseif strcmp(dose,'0dot04');d=5;
        elseif strcmp(dose,'0dot07');d=6;
        elseif strcmp(dose,'2dot40');d=7;
        end

    timepoint = M(i).timepoint;
    time = str2double(timepoint(1));
%     figure(1)
    if strcmp(species,SpeciesToPlot)
    displayname = strcat(M(i).species,'-',dose,'-',timepoint,'-',M(i).pvalue);
    % scatter(time.*M(i).ones,M(i).mrna,'MarkerFaceColor','none','MarkerEdgeColor',COLORZ{d},'DisplayName',displayname);hold on
%     scatter(time.*M(i).ones,M(i).mrna,'MarkerFaceColor','none','MarkerEdgeColor',COLORZ{d},'DisplayName',displayname);hold on
    end
end


%%%%%%%%%%%%%%%%%%%%%%%
%plotSpread
figure(96)
doses = vertcat({M.dose});
doselist = unique(doses);
specie = vertcat({M.species});
specielist = unique(specie);
times = vertcat({M.timepoint});
timelist = unique(times);
datos = vertcat({M.date});
datelist = unique(datos);

%set the bin sizes
specieidx = strcmp(specie,SpeciesToPlot);
mesRNA = vertcat({M.mrna});
allmrna = vertcat(mesRNA{specieidx});
    binN = 50;
    lowperc = 0.02;
    highperc = 0.98;
    BinCenters = linspace(prctile(allmrna,lowperc.*100),prctile(allmrna,highperc.*100),binN);

cmap = colormap(lines(50));

%correction step that shouldn't be necessary
dosss = strcmp(doses,'0dot00');
times(dosss) = {'0hr'};
%

mspread=[];
%plot for each date  (overlapping doses)
for kk = 1:length(datelist)
    Dato = datelist{kk};
for k = 1:length(specielist)
    Specie = specielist{k};
    

    cycle=1;
    for j = 1:length(timelist(1:2))%just 0hr and 1hr
        for i = 1:length(doselist)
                Dose = doselist{i};
                doseidx = strcmp(doses,Dose);
    %             Specie = SpeciesToPlot;
                specieidx = strcmp(specie,Specie);
                Time = timelist{j};
                timeidx = strcmp(times,Time);
                datoidx = strcmp(datos,Dato);

                sdidx = specieidx & doseidx & timeidx & datoidx;

                msRNA = mesRNA(sdidx);
                mRNA = vertcat(msRNA{:});
                if ~isempty(mRNA)
                    mspread{cycle} = mRNA;

                    o = M(i).ones;
                    s = o*M(i).species;
                    ss = char(s);
                    d = o*M(i).dose;
                    dd = char(d);
                    t = o*M(i).timepoint;
                    tt = char(t);
                    m = M(i).mrna;

                    std = horzcat(ss,dd,tt);

%                     cmapz{cycle} = cmap(i,:);
                    cmapz{cycle} = cmap(kk,:);
                    xvalues(cycle) = i;
                    xnames{cycle} = strcat(Dose,'-',Specie,'-',Time,'-',Dato);

    %                 plotSpread(mspread{cycle},'showMM',0,'XValues',cycle, ... 
    %                     'distributionColors',cmapz(cycle,:),'categoryMarkers',{'s'},...
    %                     'xNames',{strcat(Dose,'-',Specie,'-',Time)});
    %                 

    %                 distributionPlot(mspread{cycle},'histOpt',0,'divFactor',BinCenters,'globalNorm',3,'showMM',0,'histOri','right','XValues',cycle,'color',cmapz(cycle,:));
                    cycle=cycle+1;
                end

        end
    end

if ~isempty(mspread)
figure(33)
subplot(3,ceil(length(specielist)./2),k);hold on
plotSpread(mspread,'showMM',6,'XValues',xvalues, ... 
                    'distributionColors',cmapz,'categoryMarkers',{'s'},...
                    'xNames',xnames);
h=gca;
h.XTickLabelRotation = 90;
end

mspread=[];
cmapz=[];
xvalues=[];
xname=[];


stophere=1;
% figure(44)
% subplot(1,length(specielist),k);
% distributionPlot(mspread,'histOpt',0,'divFactor',BinCenters,'globalNorm',3,'showMM',0,'histOri','right','color',cmapz,'xNames',xnames,'xValues',xvalues)
end
end

%plot for each date (all dose and time pairs)
for k = 1:length(specielist)
    Specie = specielist{k};
    
    cycle=1;
    for j = 1:length(timelist)
        for i = 1:length(doselist)
                Dose = doselist{i};
                doseidx = strcmp(doses,Dose);
    %             Specie = SpeciesToPlot;
                specieidx = strcmp(specie,Specie);
                Time = timelist{j};
                timeidx = strcmp(times,Time);
                                
                sdidx = specieidx & doseidx & timeidx;

                msRNA = mesRNA(sdidx);
                mRNA = vertcat(msRNA{:});
                if ~isempty(mRNA)
                    mspread{cycle} = mRNA;

                    o = M(i).ones;
                    s = o*M(i).species;
                    ss = char(s);
                    d = o*M(i).dose;
                    dd = char(d);
                    t = o*M(i).timepoint;
                    tt = char(t);
                    m = M(i).mrna;

                    std = horzcat(ss,dd,tt);

                    cmapz{cycle} = cmap(i,:);
%                     cmapz{cycle} = cmap(kk,:);
                    xvalues(cycle) = cycle;
                    xnames{cycle} = strcat(Dose,'-',Specie,'-',Time);

    %                 plotSpread(mspread{cycle},'showMM',0,'XValues',cycle, ... 
    %                     'distributionColors',cmapz(cycle,:),'categoryMarkers',{'s'},...
    %                     'xNames',{strcat(Dose,'-',Specie,'-',Time)});
    %                 

    %                 distributionPlot(mspread{cycle},'histOpt',0,'divFactor',BinCenters,'globalNorm',3,'showMM',0,'histOri','right','XValues',cycle,'color',cmapz(cycle,:));
                    cycle=cycle+1;
                end

        end
    end
    
    figure(335)
subplot(2,length(specielist),k);
plotSpread(mspread,'showMM',4,'XValues',xvalues, ... 
                    'distributionColors',cmapz,'categoryMarkers',{'s'},...
                    'xNames',xnames);
h=gca;
h.XTickLabelRotation = 60;

if ~isempty(mRNA)
clear mspread cmapz xvalues xname
end

stophere=1;
% figure(44)
% subplot(1,length(specielist),k);
% distributionPlot(mspread,'histOpt',0,'divFactor',BinCenters,'globalNorm',3,'showMM',0,'histOri','right','color',cmapz,'xNames',xnames,'xValues',xvalues)
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this creates a plot of a correlation between two mrna species
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pvalues = vertcat({M.pvalue});
[~,~,~,d] = regexp(specielist,'(snail|pai1|smad7|pmepa1|tieg|bhlhe40|ctgf|wnt9a)');
uspecielist = unique(cellfun(@(x) x{1},d,'UniformOutput',0)); %this should be only ctgf and wnt9a
uspecielist = SpeciesToCorrelate;

mspec=struct();
for uspstr = uspecielist
    usp = char(uspstr);
    mspec.(usp) = [];
end

for uspstr = uspecielist
    usp = char(uspstr);
    [~,~,~,d] = regexp(specie,usp);
    uspecieidx = cellfun(@isempty,d,'UniformOutput',1);
    for datostr = datelist
        dato = char(datostr);
        datidx = strcmpi(datos,dato);
        pvs = pvalues(datidx);
        pvlist = unique(pvs);
        for pvstr = pvlist
            pv = char(pvstr);
            pvidx = strcmpi(pvalues,pv);
            mforeach = mspec.(usp);
            disp(length(find(datidx & pvidx & uspecieidx)==1))
            mspec.(usp)= vertcat(mforeach,M(datidx & pvidx & uspecieidx).mrna);
        end
    end
end
figure(5432)
fnames = fieldnames(mspec);
scatter(mspec.(fnames{1}),mspec.(fnames{2}));
[rrr,~] = corr((mspec.(fnames{1})),(mspec.(fnames{2})),'type','Pearson','rows','all');
title(strcat('correlation between-',fnames{1},'-and-',fnames{2},'-is-', num2str(rrr)));
% plotSpread(mspread);
%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%make a badass box plot!
figure(44)

for i = 1:length(M)
    o = M(i).ones;
    s = o*M(i).species;
    ss = char(s(:,1:7));
    d = o*M(i).dose;
    dd = char(d);
    t = o*M(i).timepoint;
    tt = char(t);
    m = M(i).mrna;

    std = horzcat(ss,dd,tt);

    if i==1
        mc = m;
        stdc = std;
    else


        if ~(size(std,2) == size(stdc,2))
        std(:,end+1) = char(ones(size(std(:,end))));
        end


    stdc = vertcat(stdc,std);
    mc = vertcat(mc,m);
    end
end


did = {M.dose};
tid = {M.timepoint};
dossidx = strcmp(did,'0dot00');
tid(dossidx) = {'0hr'};
for i=1:length(M)
   M(i).timepoint = tid{i}; 
end


boxplot(mc,stdc);
h=gca;
h.XTickLabelRotation=45;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%concatenate timepoints and doses
DOSESlist = vertcat({M.dose});
TIMEPOINTSlist = vertcat({M.timepoint});
SPECIESlist = vertcat({M.species});

DOSESpresent = findNumberOfVarsInSequence(DOSESlist,'[0-9]dot[0-9]+');
TIMEPOINTSpresent = findNumberOfVarsInSequence(TIMEPOINTSlist,'[0-9]hr');
SPECIESpresent = findNumberOfVarsInSequence(SPECIESlist,'(594_snail|594_smad7|594_pai1|647_snail|647_wnt9a|647_ctgf|594_wnt9a|594_ctgf|647_smad7|647_pai1|594_bhlhe40)');
Mcat=struct();
q=1;
for i=1:length(DOSESpresent)
    for j=1:length(TIMEPOINTSpresent)
        for k=1:length(SPECIESpresent)
   dose = DOSESpresent{i};
   doseLogical = strcmp(DOSESlist,dose);
   
   time = TIMEPOINTSpresent{j};
   timeLogical = strcmp(TIMEPOINTSlist,time);
   
   species = SPECIESpresent{k};
   speciesLogical = strcmp(SPECIESlist,species);
   
   onetoend = 1:length(timeLogical);
   positionsToUse = onetoend(doseLogical & timeLogical & speciesLogical);
   
   squishStructure = vertcat(M(positionsToUse));
   mrnacat = vertcat(squishStructure.mrna);
   pcat = horzcat(squishStructure.pvalue);
   
   if isempty(mrnacat)
       stophere=1;
   end
   
   if ~isempty(mrnacat)
    Mcat(q).time = time;
    Mcat(q).species = species;
    Mcat(q).dose = dose;
    Mcat(q).mrna = mrnacat;
    Mcat(q).pvalue = pcat;
    q=q+1;
   end
        end
    end
end



%%%%plot median mrna
stophere=1;
histStruct=struct();
hisMat = 1;

for q=1:length(Mcat)
   timepoint =Mcat(q).time;
   time = str2double(timepoint(1));
   
   species = Mcat(q).species;
   dose = Mcat(q).dose;
   mrna = Mcat(q).mrna;
   pvalues = Mcat(q).pvalue;
    
    medianmrna = nanmedian(mrna);
    meanmrna = nanmean(mrna);
    stdmrna = nanstd(mrna);
    cvmrna = stdmrna./meanmrna;
    sortmrna = sort(mrna);
    top = 90;
    bottom = 10;
    sortTOP = sortmrna(round(length(sortmrna).*top./100));
    sortBOT = sortmrna(round(length(sortmrna).*bottom./100));
    
    
    if strcmp(dose,'0dot00');d=1;
elseif strcmp(dose,'0dot01');d=2;
elseif strcmp(dose,'0dot02');d=3;    
elseif strcmp(dose,'0dot03');d=4; 
elseif strcmp(dose,'0dot04');d=5;
elseif strcmp(dose,'0dot07');d=6;
elseif strcmp(dose,'2dot40');d=7;
    end

    figure(2)
%     displayname = strcat(species,'-',timepoint,'-',dose,'-',pvalues);
%     errorbar(time,medianmrna,medianmrna-sortBOT,medianmrna+sortTOP,'o','MarkerEdgeColor',COLORZ{d},'DisplayName',displayname);hold on
%     title(strcat('median...',num2str(top),'%-',num2str(bottom),'% CI'));
%     xlim([-1 7])

    displayname = strcat(species,'-',timepoint,'-',dose,'-',pvalues);
    errorbar(time,medianmrna,stdmrna,'o','MarkerEdgeColor',COLORZ{d},'DisplayName',displayname);hold on
    title(strcat('median...+/- 1 std'));
    xlim([-1 7])
    
% 
%     figure(3)
%     BinEdges = 0:10:500;
%     [N,Edges,Bin] = histcounts(mrna,BinEdges);
%     barEdges = ((Edges(2)-Edges(1))/2):(Edges(2)-Edges(1)):(Edges(end)-(Edges(end)-Edges(end-1))/2);
%     
%     
%     
    
end

end




function M = runMRNA(cfile,Date)

M=struct();
M(length(cfile)).date = [];
M(length(cfile)).species = [];
M(length(cfile)).dose =[];
M(length(cfile)).timepoint = [];
M(length(cfile)).pvalue = [];
M(length(cfile)).mrna = [];
M(length(cfile)).ones = [];




parfor i=1:length(cfile)
    filename = char(cfile{i});
    [aa,bb] = regexp(filename,'(594_snail|594_smad7|594_pai1|647_snail|647_wnt9a|647_ctgf|594_wnt9a|594_ctgf|647_smad7|647_pai1|594_pmepa1|594_tieg|594_bhlhe40)');
    M(i).species = filename(aa:bb);
    % [aa,bb] = regexp(filename,'[0-9]dot[0-9]+');
    % M(i).dose = filename(aa:bb);
    [aa,bb] = regexp(filename,'(0dot00|0dot01|0dot02|0dot03|0dot04|0dot07|2dot40|low|med|high)');
    do = filename(aa:bb);


        if strcmp(do,'med')
            D = '0dot07';
        elseif strcmp(do,'high')
            D = '2dot40';
        else
            D = do;
        end
%     disp(D);
    M(i).dose = D;
    [aa,bb] = regexp(filename,'[0-9]hr');
    M(i).timepoint = filename(aa:bb);
    [aa,bb] = regexp(filename,'p[0-9]+');
    M(i).pvalue = filename(aa:bb);

    M(i).date = Date(1:10);


    A = load(filename);
    CCcells = A.CCcells;
    species = M(i).species;
    disp(cfile{i})
    fieldname = strcat('mRNA',species(1:3),'inCell');
    mrna = CCcells.(fieldname);
    M(i).mrna = mrna;
    M(i).ones = ones(size(mrna));

end
stophere=1;
end


function HOURS = findNumberOfVarsInSequence(Sequence, stringzy)
jjj=1;
HOURS=[];
for cfile = Sequence
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