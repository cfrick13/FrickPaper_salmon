function plotCORRELATIONSfromStructureBOXPLOT
close all
CORR = struct();R = struct();P = struct();Rf=struct();Rn=struct();
% drawnow

iter=100;
jacknifepercent = 0.9;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% CHOOSE PLOTTING DETAILS %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% GeneChoice = '594_snail';
% GeneChoice = '647_snail';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CountType = 'mrna_all_snail_snail';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

GeneChoice = '594_smad7';
% GeneChoice = '594_bhlhe40';
% GeneChoice = '594_pai1';
% GeneChoice = '594_ctgf';
% GeneChoice = '647_ctgf';
% GeneChoice = '594_wnt9a';
% GeneChoice = '647_wnt9a';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CountType = 'mrna_all';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CountType = 'snailcounts'; %all counts using kernel5\ no longer all counts
% 
% GeneChoice = '594_bhlhe40';
% % GeneChoice = '647_snail';
% CountType = 'm2016-25p2'; %jan25-2016 only 




% CountType = 'snailcounts'; %all counts using kernel5\ no longer all counts

% GeneChoice = '594_smad7';
% CountType = 'CountsMeanPFocuzk5all'; %all counts using kernel5\ no longer all counts
% 
% GeneChoice = '594_pai1';
% CountType = 'CountsMeanPFocuzyo5'; %only new pai1 counts using kernel5
% 
% 
% GeneChoice = '594_pai1';
% CountType = 'mcount'; %only new pai1 counts (curated) using kernel5

% GeneChoice = '594_bhlhe40';
% CountType = 'm1519new'; %

% GeneChoice = '594_bhlhe40';
% CountType = 'm1519newcur'; %Dec19 curated only

% GeneChoice = '594_bhlhe40';
% % CountType = 'm1519jan12cur';%Dec15 and 19 curated (as of jan 12)
% CountType = 'm1519full';%Dec15 and 19 curated (as of jan 12)
% CountType = 'm1519fullnocfp';
% CountType = 'm2016-25'; %jan25-2016

% GeneChoice = '594_bhlhe40';
% CountType = 'm1519new4'; %

% GeneChoice = '594_pai1';
% CountType = 'mcountdense'; %only new pai1 counts (curated) using kernel5  curated to remove low density positions

% GeneChoice = '647_pai1';
% CountType = 'CountsMeanPFocusNew'; %all data counts using kernel7

% GeneChoice = '647_ctgf';
% GeneChoice = '594_wnt9a';
% CountType = 'onlym2016-02-19'; %feb19-2016 more curation

% CountType = 'mlatestcuratedFOCUS'; %feb19-2016 more curation'mlatestcurated.mat'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

correlationType = 'spearman';
LogVSLin = 'log';
outliers = 'on';
interQuartileDistanceRange = 1.5;

msize=10;
DATAdir = 'D:\Users\zeiss\Documents\MATLAB\AllImagingDataCompiled';
cd(DATAdir)
filenamez = strcat('GeneScalarStructure',GeneChoice,CountType,'.mat');
load(filenamez); %this loads ScalarStruct, NAMO and COLOURS.
% cd('\Users\frick\Documents\MATLAB\')
ADir = 'D:\Users\zeiss\Documents\';
cd(strcat(ADir,'MATLAB'));

COLOUR = {[0.1 0.1 0], [0.2 0.2 0.0], [0.3 0.3 0.1], [0.1 0.2 0.4], [0.1 0.3 0.3], [0.0 0.2 0.4], [0.0 0.2 0.2], [0.0 0.0 0.4]}; 
lowcolor = [40  45 48]./255;
medcolor = [220 150 0]./255;
maxcolor = [255 50 0]./255;


figone = 13;

ScalarStruct = GeneScalarStruct.(strcat('i',GeneChoice,'i'));
NAMO{end+1} = 'mRNA';
DoseNamess = fieldnames(ScalarStruct);
DoseNamessort = sort(DoseNamess');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% determine which doses to remove %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% removedim = 1;
% removedim = length(DoseNamessort);
removedim = [];
    cycledim = 1:length(DoseNamessort);
    cycledim(removedim)=[];
    DoseNames = cell(1,length(cycledim));
    cycle=1;
        for i = cycledim
        DoseNames{cycle} = DoseNamessort{i};
        cycle=cycle+1;
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    SCALAR = cell(1,length(DoseNames));
    SCALARsansOL = cell(1,length(DoseNames));
    DoseStr = cell(1,length(DoseNames));
    i=1;
    for Doses = sort(DoseNames)
        dosestring = Doses{1};
            
            %determines color of points to be plotted.
            if str2double(dosestring(end-2:end-1))<3
                COLOUR{i} = lowcolor;
            elseif str2double(dosestring(end-2:end-1))<7
                COLOUR{i} = medcolor;
            else
                COLOUR{i} = maxcolor;
            end

        SCALAR{i} = ScalarStruct.(char(Doses));
        scalarcounts = ScalarStruct.(char(Doses));
        emmiescount = scalarcounts(1,~isnan(scalarcounts(end,:)));
        DoseStr{i} = char(Doses);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% REMOVE OUTLIERS (determine idx of outliers)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if strcmp(outliers,'off')
        fulloutlierLogical = zeros(size(scalarcounts(end,:)));
        else
            shimrnasub = scalarcounts(end,:);
            if strcmp(LogVSLin,'log')
            [~,~,fulloutlierLogical] = removeOutliers(log10(shimrnasub),[],interQuartileDistanceRange);
            else
            [~,~,fulloutlierLogical] = removeOutliers(shimrnasub,[],interQuartileDistanceRange);
            end
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        shimmy = ScalarStruct.(char(Doses));
        if i==1
            ScalarConcatenate = shimmy(:,~fulloutlierLogical);
        else
            ScalarConcatenate = horzcat(ScalarConcatenate,shimmy(:,~fulloutlierLogical));
        end
        SCALARsansOL{i} = shimmy(:,~fulloutlierLogical);
                
            i=i+1;
            disp(strcat('Cells with Counts=',char(Doses),'...',num2str(size(emmiescount,2))));
    end

    
    
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Determine Correlations %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ScalarConcatenateNANremoved = ScalarConcatenate(:,~isnan(ScalarConcatenate(end,:)));

%remove outliers from the ScalarConcatenated Dataset
Shimcat = ScalarConcatenateNANremoved(1:end-1,:);
shimrna = ScalarConcatenateNANremoved(end,:);

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% %%% isolate the dynamic region of data
%         fchimcat = Shimcat(5,:);
%         midx = shimrna>30 & shimrna<474;
%         fcidx = fchimcat<3.5; % 
%         idx = midx&fcidx;
%         Shimcat = Shimcat(:,idx);
%         shimrna = shimrna(idx);
%         disp('dynamic isolation')
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    
    rr=zeros(size(Shimcat,1),iter);
    pp=zeros(size(Shimcat,1),iter);
    rrK=zeros(size(Shimcat,1),iter);
    ppK=zeros(size(Shimcat,1),iter);
    rrS=zeros(size(Shimcat,1),iter);
    ppS=zeros(size(Shimcat,1),iter);
    rrmse=zeros(size(Shimcat,1),iter);

    
        [rrr,ppp] = corr(Shimcat',shimrna','type','Pearson','rows','all');
        [rrrK,pppK] = corr(Shimcat',shimrna','type','Kendall','rows','all');
        [rrrS,pppS] = corr(Shimcat',shimrna','type','Spearman','rows','all');
  
    parfor i = 1:iter
 
        shidx = randi(size(Shimcat,2),1,round(size(Shimcat,2).*jacknifepercent));
        Shimcatsub = Shimcat(:,shidx);
        shimrnasub = shimrna(:,shidx);

        
        [r,p] = corr(Shimcatsub',shimrnasub','type','Pearson','rows','all');
        [rK,pK] = corr(Shimcatsub',shimrnasub','type','Kendall','rows','all');
        [rS,pS] = corr(Shimcatsub',shimrnasub','type','Spearman','rows','all');

        mdlr = zeros(1,size(Shimcatsub,1));
        mdlp = zeros(1,size(Shimcatsub,1));
        rmse = zeros(1,size(Shimcatsub,1));
        
%             for y = 1:size(Shimcatsub,1)
%                 mdl = fitlm(Shimcatsub(y,:),shimrnasub);
%                 mdlr(y) = sqrt(mdl.Rsquared.Adjusted); 
%                 mdlp(y) = mdl.Coefficients.pValue(2);
%                 rmse(y) = mdl.RMSE;
%             end
%          
            disp(i)
            
        rr(:,i) = r;
%         rr(:,i) = r;
        pp(:,i) = p;
        rrK(:,i) = rK;
        ppK(:,i) = pK;
        rrS(:,i) = rS;
        ppS(:,i) = pS;
        rrmse(:,i) = rmse;
        
    end
    
%     rstderror = std(rr,[],2)./sqrt(sum((~isnan(Shimcat)),2));
    
        r = nanmean(rr,2);%%%%%
        rstd = nanstd(rr,[],2);
        rf = prctile(rr,5,2);
        rn = prctile(rr,95,2);
        rstderror = rstd./sqrt(round(size(Shimcat,2).*jacknifepercent));
    rK = nanmean(rrK,2);
    rKstd = nanstd(rrK,[],2);
    rKf = prctile(rrK,5,2);
    rKn = prctile(rrK,95,2);
        rS = nanmean(rrS,2);
        rSstd = nanstd(rrS,[],2);
        rSf = prctile(rrS,5,2);
        rSn = prctile(rrS,95,2);
    p = nanmean(pp,2);
    pK = nanmean(ppK,2);
    pS = nanmean(ppS,2);
        
        
    rmse = nanmean(rrmse,2);
    rstdermse = nanstd(rrmse,[],2)./sqrt(round(size(Shimcat,2).*jacknifepercent));
    
    

% R.pearson=r;
R.pearson=rrr;
R.kendall=rrrK;
R.spearman=rrrS;

Rf.pearson=rf;
Rf.pearson=rf;
Rf.kendall =rKf;
Rf.spearman = rSf;

Rn.pearson=rn;
Rn.kendall =rKn;
Rn.spearman = rSn;
% R.mdlr = mdlr;

P.pearson=ppp;
P.kendall=pppK;
P.spearman=pppS;
% P.mdlr = mdlp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






[a,b] = regexp(NAMO,'[0-9]+');
    for i=1:length(a)
        aa = a{i};
        bb = b{i};
        if ~isempty(aa)
            namoidx = aa:bb;
            namo = NAMO{i};
            namo(namoidx)=[];
            NAMOsansNumbers{i} = namo;
        else
            NAMOsansNumbers{i} = NAMO{i};
        end
        
        while ~isempty(regexp(NAMOsansNumbers{i},'\-*'));
        [c,d] = regexp(NAMOsansNumbers{i},'\-*');
        namo = NAMOsansNumbers{i};
        namo(c:d) = [];
        NAMOsansNumbers{i} = namo;
        end
    end
    
uniqueNAMO = unique(NAMOsansNumbers,'stable');
% 23:36 are all of the unique NAMO that have multiple time points

    cycle=1;
    for i = 23:length(uniqueNAMO)-1
        uniqueVarList{cycle} = find(strcmp(NAMOsansNumbers,uniqueNAMO{i})==1);
        cycle=cycle+1;
    end


cy=1;
for b=uniqueVarList;
    a=b{1};
    CORR(cy).pearson = find(r(a)==max(r(a)),1,'last')+a(1)-1;
    CORR(cy).spearman = find(rS(a)==max(rS(a)),1,'last')+a(1)-1;
    CORR(cy).kendall = find(rK(a)==max(rK(a)),1,'last')+a(1)-1;
%     CORR(cy).mdlr = find(mdlr(a)==max(mdlr(a)),1,'last')+a(1)-1;
    
    CORR(cy).pearsonP = find(p(a)==min(p(a)),1,'last')+a(1)-1;
    CORR(cy).spearmanP = find(pS(a)==min(pS(a)),1,'last')+a(1)-1;
    CORR(cy).kendallP = find(pK(a)==min(pK(a)),1,'last')+a(1)-1;
%     CORR(cy).mdlrP = find(mdlp(a)==min(mdlp(a)),1,'last')+a(1)-1;
    cy=cy+1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%  BOXPLOT  %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=1;
b=0;
c=1;
NameCell = cell(1,size(ScalarConcatenateNANremoved,2));
mCat = nan(1,size(ScalarConcatenateNANremoved,2));
    for Doses = sort(DoseNames)

        %%%%%%%%%%
        mum = ScalarStruct.(char(Doses));
        mums = mum(end,:);
        mumsie = mums(~isnan(mums));
        MMlength = length(mumsie);

        %%%%%%%%%%
        shimrnasub = mumsie;
        if strcmp(LogVSLin,'log')
            [~,~,fulloutlierLogical] = removeOutliers(log10(shimrnasub),[],interQuartileDistanceRange);
        else
            [~,~,fulloutlierLogical] = removeOutliers(shimrnasub,[],interQuartileDistanceRange);
        end

        %%%%%%%%%%
        for i = a:a+MMlength-1
            NameCell{b+i} = char(Doses);
            if strcmp(LogVSLin,'log')
                mCat(b+i) = log10(mumsie(i));
            else
                mCat(b+i) = mumsie(i);    
            end
        end

        %%%%%%%%%%
        b=MMlength+b;
        if strcmp(LogVSLin,'log')
            MSM{c}=log10(mumsie);
            MSMout{c} = log10(mumsie(fulloutlierLogical));
        else
            MSM{c}=mumsie;
            MSMout{c} = mumsie(fulloutlierLogical);
        end
        
        %%%%%%%%%%
        if ~isempty(mumsie)
            c=c+1;
        end
    end

    figure(789)
    if ~(sum(cellfun(@isempty,MSMout))==length(MSMout))
        plotSpread(MSMout);hold all
        h=gca;
            for m = h.Children'
                m.Color = 'r';
                m.Marker = 'o';
                m.MarkerSize = 12;
            end
    end

plotSpread(MSM);hold all
boxplot(mCat,NameCell);hold all
h=gca;

    if strcmp(LogVSLin,'log')
    h.YTickLabel = 10.^h.YTick;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









%%%%%
featuresToPlot = horzcat(CORR.(correlationType));
choosefeaturefortime = 5;
idx = find(uniqueVarList{choosefeaturefortime} == featuresToPlot(choosefeaturefortime));
for i=1:length(featuresToPlot)
   varlist = uniqueVarList{i};
   if idx>length(varlist)
       idxi = length(varlist);
   else
       idxi=idx;
   end
    featuresToPlot(i) =  varlist(idxi);
end






%%%%%%
%MAKE time correlation errorbar PLOT of time correlations!!!!
%%%%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% DETERMINE COLOR MAP FOR PLOTTING TRACES %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:length(r)
    
    if ~isempty(regexp(NAMO{i},'\-B'));
%         cmap{i} = 'y';;
        cmap{i} = [0.6 0.6 0]; 
    elseif ~isempty(regexp(NAMO{i},'\/B'));
        cmap{i} = [0.2 0.6 0] ;
    elseif ~isempty(regexp(NAMO{i},'\Time*'));
        cmap{i} = [0 0.6 0.6];
    elseif ~isempty(regexp(NAMO{i},'\Percent*'));
        cmap{i} = [0.6 0.6 0.6];
    else
        cmap{i} = [0.6 0.1 0.0];
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



figure(9182)   
cpr = R.(correlationType);
cprf = cpr-Rf.(correlationType);
cprn = Rn.(correlationType)-cpr;
% cpr = r;
% cpr = rmse;
cyc=1;
    for i=1:length(uniqueVarList)
       varlist = uniqueVarList{i};
       if sum(factor(i-1)==3)>0
           cyc=cyc+1;
       end
       subplot(2,3,cyc);
       plot(cpr(varlist),'Color',cmap{varlist(1)},'LineWidth',2,'Marker','s','MarkerSize',6,'MarkerEdgeColor','none','MarkerFaceColor',cmap{varlist(1)});hold on
%        e = errorbar(cpr(varlist),rstd(varlist),'Color','k','LineWidth',1);hold on
       e = errorbar(cpr(varlist),rstderror(varlist),'Color','k','LineStyle','none','LineWidth',1);hold on
%          e = errorbar([1:16],cpr(varlist),rf(varlist),rn(varlist),'Color','k','LineStyle','none','LineWidth',1);hold on
% %        e.DisplayName = uniqueNAMO{i+21};
    
h=gca;
h.YLabel.String = 'Coefficient of Determination, R^2';
h.XLabel.String = 'Minutes';
h.XTick = [1 6 11 16];
h.XTickLabel = [0 20 40 60];
h.XLim = [0 16];
h.YLim = [-0.02 0.3];
h.Title.String = NAMOsansNumbers{varlist(1)};
% axis('tight')
h.Box = 'off';
    end
   
    
    
figure(9183)   
cpr = R.(correlationType);

% cpr = r;
% cpr = rmse;
cyc=1;
    for i=1:length(uniqueVarList)
       varlist = uniqueVarList{i};
       if sum(factor(i-1)==3)>0
           cyc=cyc+1;
       end
       subplot(2,3,cyc);
       plot(cpr(varlist),'Color',cmap{varlist(1)},'LineWidth',2,'Marker','s','MarkerSize',6,'MarkerEdgeColor','none','MarkerFaceColor',cmap{varlist(1)});hold on
%        e = errorbar(cpr(varlist),rstd(varlist),'Color','k','LineWidth',1);hold on
       e = errorbar(cpr(varlist),rstdermse(varlist),'Color','k','LineStyle','none','LineWidth',1);hold on
%        e.DisplayName = uniqueNAMO{i+21};
    
h=gca;
h.YLabel.String = 'RMSE';
h.XLabel.String = 'Minutes';
h.XTick = [1 6 11 16];
h.XTickLabel = [0 20 40 60];
h.XLim = [0 16];
h.YLim = [50 150];
h.Title.String = NAMOsansNumbers{varlist(1)};
% axis('tight')
h.Box = 'off';
   end








%%%%%%%%%%%%%%%%%%
%Make  awesome barplot of correlations (no time variables)
%%%%%%%%%%%%%%%%%%

[numpresent,~] = regexp(NAMO,'[0-9]');
nump = cellfun(@isempty,numpresent);
nmp = find(nump==0,1,'first');

for i = 1:nmp-1
    
    if ~isempty(regexp(NAMO{i},'\-B'));
%         cmap{i} = 'y';;
        cmap{i} = [0.6 0.6 0]; 
    elseif ~isempty(regexp(NAMO{i},'\- B'));
        cmap{i} = [0.6 0.6 0]; 
    elseif ~isempty(regexp(NAMO{i},'\/B'));
        cmap{i} = [0.2 0.6 0] ;
    elseif ~isempty(regexp(NAMO{i},'\Time*'));
        cmap{i} = [0 0.6 0.6];
    elseif ~isempty(regexp(NAMO{i},'\Percent*'));
        cmap{i} = [0.6 0 0.6];
    else
        cmap{i} = [0.4 0.4 0.4];
    end
    
end


figure(91899)
for i = 1:nmp-1
   bar(i,cpr(i),'FaceColor',cmap{i});hold on
end
%    errorbar(r(1:nmp-1),rstd(1:nmp-1),'LineStyle','none','Color','k')
   errorbar(1:length(cpr(1:nmp-1)),cpr(1:nmp-1),cprf(1:nmp-1),cprn(1:nmp-1),'LineStyle','none','Color','k')
h=gca;
h.XTick = 1:nmp-1;
h.XTickLabel = NAMO(1:nmp-1);
h.XTickLabelRotation = 60;
h.YLabel.String = 'Coefficient of Determination, R^2';
h.YLim = [0 0.7];   
h.XLim = [0 nmp];

    
stophere=1;






%%%%%


%%%%%%%%%% FIGURE OF ALL THE CORRELAITONS
% % % figure(figthr)
% % % coloror = vertcat(COLORS');
% % % colorors=cell2mat(coloror);
% % % area = ones(1,length(r)).*4;
% % % area(horzcat(CORR.pearson)) = 15;
% % % subplot(2,3,1);
% % %     scatter(1:length(r),r,area,colorors);
% % %     title('Pearson corrcoeff')
% % %     ylim([0 1])
% % % subplot(2,3,4);
% % %     scatter(1:length(p),p,area,colorors);
% % %     title('Pearson pvalue')
% % %     ylim([0 median(p)*2])
% % % subplot(2,3,2);
% % %    scatter(1:length(rK),rK,area,colorors);
% % %     title('Kendall corrcoeff')
% % %     ylim([0 1])
% % % subplot(2,3,5);
% % %     scatter(1:length(pK),pK,area,colorors);
% % %     title('Kendall pvalue')
% % %     ylim([0 median(pK).*2])
% % % subplot(2,3,3);
% % %     scatter(1:length(rS),rS,area,colorors);  
% % %     title('Spearman corrcoeff')
% % %     ylim([0 1])
% % % subplot(2,3,6);
% % %     scatter(1:length(pS),pS,area,colorors);
% % %     title('Spearman pvalue')
% % %     ylim([0 median(pS).*2])
%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%PLOT CORRELATIONS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
jj=1;
% featuresToPlot = [32 46 60];
% featuresToPlot = [1 2 3 4 5];

% featuresToPlot = [32,46,60,74,88,102,116,130,144,158,170,182,230,240,250,260,272,284];
minutes =1;
featuresToPlot = [32,46,60,72,88,100,116,130,144,158,170,182,230,240,250,260,272,284] +minutes;
% featuresToPlot = [32 46 60];
    for chooseScalar = featuresToPlot
    
    Nvalue=0;
        for i=1:length(SCALAR)
            scalarSub = SCALAR{i};

            figure(figone)
            ax = subplot(4,round(length(featuresToPlot)./4)+1,jj);

            shimrnasub = scalarSub(end,:);
            shimcat = scalarSub(chooseScalar,:);
                if strcmp(LogVSLin,'log')
                    [~,~,fulloutlierLogical] = removeOutliers(log10(shimrnasub),log10(shimcat),interQuartileDistanceRange);
                    shimrnaout = shimrnasub(~fulloutlierLogical);
                    shimcatout = shimcat(~fulloutlierLogical);
                else
                    [shimrnaout,shimcatout,~] = removeOutliers(shimrnasub,shimcat,interQuartileDistanceRange);
                end
                    if strcmp(outliers,'off')
                    h = plot(shimcat,shimrnasub);hold on
                    Nvalue = length(shimrnasub(~isnan(shimrnasub)))+Nvalue;
                    else
                        if ~isempty(shimrnaout)
                            shimidx = ~(shimcatout>0); shimcatout(shimidx)=[];shimrnaout(shimidx)=[];
                            h = plot(shimcatout,shimrnaout);hold on
                            h.DisplayName = DoseStr{i};
                        end
                    Nvalue = length(shimrnaout(~isnan(shimrnaout)))+Nvalue;
                    end
         
            h.LineStyle = 'none';
            h.Marker = 'o';
            h.MarkerSize = msize;
            h.MarkerFaceColor = 'none';
            h.MarkerEdgeColor = COLOUR{i};
%             disp(DoseStr{i});
            axes(ax);
            ax.YLabel.String = 'mRNA transcripts';
            ax.XLabel.String = NAMO{chooseScalar};
            ax.YScale = LogVSLin;
        %     m.YScale = 'linear';
            ax.XScale = LogVSLin;  
        end
        rcr = R.(correlationType);
        titleSTR{1} = NAMO{chooseScalar};
            [ag,~] = regexp(GeneChoice,'\_*');
            GeneChoice(ag) = '-';
        titleSTR{2} = strcat(GeneChoice);
        titleSTR{3} = strcat('N = ', num2str(Nvalue),'     -    R = ',num2str(round(rcr(chooseScalar),2,'significant')));
    title(char(titleSTR));
    
    [sortedscal,~] = sort(ScalarConcatenateNANremoved(chooseScalar,:));
    sortedms = sort(ScalarConcatenateNANremoved(end,:));

    %%%%%%%%%%%%%
   
    Lowscal = sortedscal(ceil(length(sortedscal).*0.05));
    if Lowscal == 0
        axis('tight')
    else
 
        Lowms =  sortedms(ceil(length(sortedms).*0.02));
        Highscal = sortedscal(round(length(sortedscal).*0.95));
        Highms =  sortedms(round(length(sortedms).*0.95));
            if strcmp(LogVSLin,'log')
                if Lowscal<0
                    xtick = logspace(log10(1),log10(Highscal),5);
                    xtick = sort(horzcat(xtick,Lowscal));
                else
                    xtick = logspace(log10(Lowscal),log10(Highscal),5);
                end
                    ytick = logspace(log10(Lowms),log10(Highms),5);
            else
                xtick = linspace((Lowscal),(Highscal),5);
                ytick = linspace((Lowms),(Highms),5);
            end

            ax.XTick = round(xtick,2,'significant');
            ax.XTickLabel = round(xtick,2,'significant');
            ax.YTick = round(ytick,2,'significant');
            ax.YTickLabel = round(ytick,2,'significant');
                if Lowscal<0
                    ax.XLim = [min(xtick).*5 Highscal.*1.5];
                else
                    ax.XLim = [Lowscal./1.5 Highscal.*1.5];
                end
            ax.YLim = [Lowms Highms.*1.5];
    end
        
    jj=jj+1;
    end
    
    
figure(figone)
stuffgoodforplotting
stophere=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


function [newDataSet,newCompanionSet,fulloutlierLogical] = removeOutliers(DataSet,CompanionSet,InterQuartileDistanceRange)
% function [newDataSet,newCompanionSet] = removeOutliers(DataSet,CompanionSet,InterQuartileDistanceRange)

if isempty(InterQuartileDistanceRange)
InterQuartileDistanceRange=1.5;
end

dataset  = DataSet;
datasetnums = dataset(~isnan(dataset));

[datasetsort,sortLog] = sort(datasetnums);
quartiles = quantile(datasetsort,[0.25 0.5 0.75]);
IQD = abs(quartiles(3)-quartiles(1));
outlierLogical = (datasetsort > quartiles(3)+IQD.*InterQuartileDistanceRange | datasetsort < quartiles(1) - IQD.*InterQuartileDistanceRange);
fulloutlierLogical = (dataset > quartiles(3)+IQD.*InterQuartileDistanceRange | dataset < quartiles(1) - IQD.*InterQuartileDistanceRange);


newDataSet = datasetsort(~outlierLogical);

if isempty(CompanionSet)
    newCompanionSet = [];
else
    CompanionSetnums = CompanionSet(~isnan(dataset));
    CompanionSetsort = CompanionSetnums(sortLog);
    newCompanionSet = CompanionSetsort(~outlierLogical);
end
    
end



function stuffgoodforplotting


% cmapfinal = cmapgenerator;
% imshow(reshape(cmapfinal,[1 size(cmapfinal,1) size(cmapfinal,2)]));
% cmap=cmapfinal;

f=gcf;
f.Children;
    for h=f.Children'
        gca;
        axis('tight')
        Pos = h.Position;
        Pos(3:4) = [0.08 0.12];
        Pos(3:4) = [0.08 0.12]./1.5;
        h.Position = Pos;
        h.FontSize = 10;
        h.FontName = 'Avenir Next';
        llength = length(h.Children);
        cmap = cmapGenerator(llength);
        i=1;
            for l = h.Children'
%                 l.MarkerFaceColor=cmap(size(cmap,1)-i,:);
                l.MarkerFaceColor=cmap(i,:);
                l.MarkerEdgeColor = [0 0 0];
                l.MarkerSize = 5;
                i=i+1;
            end
    end
end

function cmapfinal = cmapgen(lenghtofdata)

oc = [1 2 1]; %orderofcmapmaking [luminance alpha  beta]
N = [50000 1000 8];

% cmap = colormap(hsv(1000000));
cmap = rand(1000000,3);
lab = im2double(applycform(cmap,makecform('srgb2lab')));
luminance = lab(:,oc(1)); %luminance of oc(1)=1;
[L,I] = sort(luminance);
yf = round(length(I)./2);
xf = yf-N(1);
labone = lab(I(xf:yf),:); %newlab has rouhgly same luminance values
cmapone = cmap(I(xf:yf),:);

alpha = labone(:,oc(2)); %alpha if oc(2) = 2;
[A,I] = sort(alpha);
yf = round(length(I)./2);
xf = yf-N(2);
labtwo = labone(I(xf:yf),:);
cmaptwo = cmapone(I(xf:yf),:);

beta = labtwo(:,oc(3)); %beta if oc(3) =3;
[B,I] = sort(beta); %sort based on green and magenta
idxs = round(linspace(1,length(I),N(3)));
labfinal  = labtwo(I(idxs),:);
cmapfinal = cmaptwo(I(idxs),:);

end



    