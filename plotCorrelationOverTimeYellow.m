cd('D:\Users\zeiss\Documents\MATLAB\')
% load('mrnapostSNAIL-FINALcurated100iter-notopdose.mat') %SNAIL no top dose
% load('mrnapostSNAIL-FINALcurated100iter.mat') %SNAIL
% load('mrnapostPAI1-FINALcurated100iter.mat') %PAI1
% load('mrnapostBHLHE40-FINALcurated100iter.mat') %BHLHE40
% load('mrnapostSNAIL-FINALcurated100iter-notopdose647.mat') %SNAIL+647
% load('mrnapostSNAIL-FINALcurated100iter-647.mat') %Snail+647 w/top dose
% % load('mrnapostWNT9A-preFINALcurated100iter-2016-02-17.mat') %wnt9a594+647
% load('mrnapostWNT9A-betterFINALcurated100iter-2016-02-17-b.mat')%wnt9a594+647 better curation
% % load('mrnapostCTGF-preFINALcurated100iter-2016-02-17.mat') %CTGF594+647
% load('mrnapostCTGF-betterFINALcurated100iter-2016-02-17-b.mat')%CTGF594+647 better curation
% load('mrnapostBHLHE40-FINALcurated100iter2016.mat')
% load('mrnapostSNAIL-betterFINALcurated100iter-ALL.mat')

% load('mrnapostBHLHE40-fullFINALcurated100iter-ALL.mat');
% load('mrnapostPAI1-fullFINALcurated100iter-ALL.mat');
% load('mrnapostSMAD7-fullFINALcurated100iter-ALL.mat');

% load('mrnapostCTGF-fullFINALcurated100iter-ALL.mat');
%     load('mrnapostCTGF-09-fullFINALcurated100iter-ALL.mat');
% load('mrnapostWNT9A-fullFINALcurated100iter-ALL.mat');
%     load('mrnapostWNT9A-09-fullFINALcurated100iter-ALL.mat');    
% load('mrnapostSNAIL-fullFINALcurated100iter-ALL.mat');

corrarray = {'pearson','spearman','kendall'};
for subnumber=1:length(corrarray)
    linestylez='-';
    correlationType = corrarray{subnumber};

% % correlationType = 'pearson';
% linestylez = '-';
% subnumber=1;
% 
% correlationType = 'spearman';
% linestylez = '-';
% subnumber=2;
% 
% correlationType = 'kendall';
% linestylez = ':';
% subnumber=3;

conf = 90;
lp = 50-conf./2;
hp = 50+conf./2;
%         r = nanmean(rr,2);
        r=rrr;
        rstd = nanstd(rr,[],2);
        rf = r-prctile(rr,lp,2);
        rn = prctile(rr,hp,2)-r;
        rstderror = rstd./sqrt(round(size(Shimcat,2).*jacknifepercent));
%     rK = nanmean(rrK,2);
    rK=rrrK;
    rKstd = nanstd(rrK,[],2);
    rKf = rK - prctile(rrK,lp,2);
    rKn = prctile(rrK,hp,2)-rK;
%         rS = nanmean(rrS,2);
        rS=rrrS;
        rSstd = nanstd(rrS,[],2);
        rSf = rS-prctile(rrS,lp,2);
        rSn = prctile(rrS,hp,2)-rS;

R.pearson=rrr;
R.kendall=rrrK;
R.spearman=rrrS;

Rf.pearson=rf;
Rf.kendall =rKf;
Rf.spearman = rSf;

Rn.pearson=rn;
Rn.kendall =rKn;
Rn.spearman = rSn;

figure(9183)
% cpr = R.pearson;
cpr = R.(correlationType);
cprn = Rn.(correlationType);
cprf = Rf.(correlationType);
% if strcmp(correlationType,'pearson')
%     cpr = sqrt(cpr);
% end
% cpr = rmse;




cyc=1;
cycle=0;
for i=1:length(uniqueVarList)
    cycle=cycle+1;
varlist = uniqueVarList{i};
        if sum(factor(i-1)==3)>0
        cyc=cyc+1;
        end
    if cyc ==2
    if 1==1 %~sum(factor(cycle+1)==3)>0    
%     subplot(2,3,cyc);
    subplot(2,3,subnumber);    
    if i==4
        aaa=0;
    elseif i==5
        aaa=0.15;
    elseif i==6
        aaa = 0.3;
    end
    plot([1:length(varlist)]+aaa,cpr(varlist),'Color',cmap{varlist(1)},'LineStyle',linestylez,'LineWidth',2,'Marker','s','MarkerSize',6,'MarkerEdgeColor','none','MarkerFaceColor',cmap{varlist(1)});hold on
    %        e = errorbar(cpr(varlist),rstd(varlist),'Color','k','LineWidth',1);hold on
%     e = errorbar([1:length(varlist)]+aaa,cpr(varlist),rstd(varlist),'Color','k','LineStyle','none','LineWidth',1);hold on
      e = errorbar([1:length(varlist)]+aaa,cpr(varlist),cprf(varlist),cprn(varlist),'Color',cmap{varlist(1)}./2,'LineStyle','none','LineWidth',1);hold on
%        e.DisplayName = uniqueNAMO{i+21};
    h=gca;
    h.YLabel.String = 'Correlation Coefficient';
    h.XLabel.String = 'Minutes';
    h.XTick = [1 6 11 16];
    h.XTickLabel = [0 20 40 60];
    h.XLim = [0 16];
    h.YLim = [-0.2 0.7];
%     h.Title.String = NAMOsansNumbers{varlist(1)};
    h.Title.String = correlationType;
    % axis('tight')
    h.Box = 'off';
    end
    end
end




% cpr = R.(correlationType);
% cprf = cpr-Rf.(correlationType);
% cprn = Rn.(correlationType)-cpr;
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
subplot(2,3,subnumber)
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
h.Title.String = correlationType;
    
stophere=1;


end
figure(9183)
tex = text(100,100,GeneChoice);
tex.Units = 'pixels';
tex.Position = [-150 -100 100];

