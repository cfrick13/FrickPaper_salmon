%DetermineMutualInformation(SigFeaturesArray,k,selectedFeatures,dimensions)
% see Estimating Mutual Information by local gaussian approximation by
% Alexander Kraskov, Harald St¨ogbauer, and Peter Grassberger (2008)
%   "SigFeaturesArray" is a cell array where each cell is a different condition and
%   the items contained within the cells are matrices with the experimental
%   measurement being made (such as a measure of signaling response) in
%   (dim=2). (dim=1) corresponds to the cell (or obseration# where the
%   measurement was made)
%     
%     the matrix within SigFeaturesArray is assemble as shown below
%  [measure1|observation1  measure2|observation1 ..... measureN|observation1;...
%          :                         :                           :          
%          :                         :                           :
%  [measure1|observationM  measure2|observationM ..... measureN|observationM]                                     :


%   "k" is the number of nearest neighbors
%   "selectedFeatures" defines the indices of the specific signaling features within "SigFeaturesArray" to be
%   quantified. Should be a vector of scalars.
%   "dimensions" is the dimensions to determine the information at (1d ([x]), 2d ([x,y]), or 3d ([x,y,z]) timepoint
%   combinations. Should be a scalar : 1, 2 or 3.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd('D:\Users\zeiss\Documents\MATLAB\AllImagingDataCompiled\')
% load('preinfo-2016-02-29.mat')
load('preinfo-2016-03-03.mat')
cd('D:\Users\zeiss\Documents\MATLAB\')
% load('/Users/frick/Documents/Goentoro_Lab/Writing/Information Paper/FIGURES/2016_02_29/Figure3/INFORMATIONz01-Mar-2016-1d.mat')

cycle=1;
fnames = sort(fieldnames(ScalarStructOnly))
for fn = fnames'
sc = ScalarStruct.(fn{1});
if size(sc,2)>50
INFOyo{cycle} = sc(1:end-1,:);
cycle=cycle+1;
end
end
selectedFeaturesoned = [[1:21] [65:78] [93:106]];
% Chris_Info_CalculationsUpdatedALLDMACnew(INFOyo, 3, selectedFeaturesoned,1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = DetermineMutualInformation(SigFeaturesArray,k,selectedFeatures,dimensions)
NumSignalConditionsInStruct = length(SigFeaturesArray);
iter = 100;
featuresVector = selectedFeatures;
assembledFeatures = chooseQFfunct(SigFeaturesArray,featuresVector,dimensions); %assembles an array of 1d ([x]), 2d ([x,y]), or 3d ([x,y,z]) timepoint combinations depedning on the dimension input
lengthOfFeatures = length(assembledFeatures);    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% chosenScalars = {[24 30 36]};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  
%define the smoothness and the means (mus) for the signaling input
%probability distributions
smthness=2000;
number_of_PofS=100; %number of probability distributions to test
lowsigma=0.0001;
highsigma=1;
InputSignals = generateInputDistributions(NumSignalConditionsInStruct,smthness,number_of_PofS,lowsigma,highsigma); %generates a continuum of unimodal, bimodal, trimodal inputs
%%%%%%%%%%%
InputSignals(1,:) = ones(1,size(InputSignals,2))./size(InputSignals,2); %set the first input to be a uniform input distribution. 
%%%%%%%%%%%


InfoMeanForAllResponses = zeros(size(InputSignals,1),lengthOfFeatures);
InfoStdForAllResponses = zeros(size(InputSignals,1),lengthOfFeatures);
for cyclenumber = 1:size(InputSignals,1)
    PofSvector = InputSignals(cyclenumber,:); %vector where each number is a probability of a give signal (condition) 
    disp(cyclenumber); %Display cycle to track progress

    %Define Data
    DataCell = cell(NumSignalConditionsInStruct,length(assembledFeatures)); %initialize
    PRofScellarray = makePRofScellarray(DataCell,SigFeaturesArray,assembledFeatures);  %and then determine PRofScellarray
    numberOfSignals = size(PRofScellarray,1);

    
    %determine the value of PSI(k) and PSI(N) and %Determine (H(Z|S))
    Digam = zeros(numberOfSignals,lengthOfFeatures); %initialize
    entropyRgS = zeros(numberOfSignals,lengthOfFeatures); %initialize
    Nvalue = zeros(numberOfSignals,lengthOfFeatures); %initialize
    SampSize = zeros(1,size(PRofScellarray,1))'; %initialize
    
    DigK = psi(k); %determine the value of digammak or PSI(k)
    for d = 1:numberOfSignals
        for iii=1:lengthOfFeatures
            ResponsesGivenS = PRofScellarray{d,iii};
            SampSize(d) = size(ResponsesGivenS,2);
            responsesGivenSvector =ResponsesGivenS';
            if sum(sum(isnan(responsesGivenSvector)))>0
                idx = reshape(isnan(ResponsesGivenS),size(ResponsesGivenS));
                [~,col,~] = find(idx==1);
                [~,dim2] = size(ResponsesGivenS);
                nline = 1:dim2;
                nline(col)=[];
                responsesGivenSvector = ResponsesGivenS(:,nline)';  
            end

            featuresSelected = assembledFeatures{iii};
            if dimensions >1 
                dmnsn = determineDimension(featuresSelected,responsesGivenSvector);
            else %else statement is for speed
                dmnsn=1;
            end

            cdz = (pi.^(dmnsn./2))./gamma(1+dmnsn./2)./2.^dmnsn; %volume of d-dimensional unit ball.

            if dimensions >1
                urgs = [];
                responsesGivenSvector = responsesGivenSvector;
            elseif featuresSelected==2 || featuresSelected==22;
                urgs = unique(responsesGivenSvector);
                responsesGivenSvector = urgs;
            else
                urgs = [];
                responsesGivenSvector = responsesGivenSvector;
            end

            [~, Eps] = knnsearch(responsesGivenSvector,responsesGivenSvector,'K', k);

            Eps = Eps(:,k).*2;%Eps is the distance x 2
            DigammaN = psi(length(responsesGivenSvector));
            Digam(d,iii) = DigammaN;
            entropyRgS(d,iii) = (dmnsn./length(responsesGivenSvector))*sum(log(Eps)) - DigK + DigammaN + log(cdz); %step1[use all cells in each distribution to estimate each h(R|S=s)
            Nvalue(d,iii) = length(responsesGivenSvector);
        end
    end
InfoVec = calculateInfo(PRofScellarray,SampSize,numberOfSignals,PofSvector,lengthOfFeatures,k,DigK,entropyRgS,iter,assembledFeatures,dimensions);

Infos = log2(exp(1)).*(-InfoVec); %Calculate Information  
Infos(Infos==-inf)=NaN;
InfoMeanForAllResponses(cyclenumber,:) = nanmean(Infos,1);
InfoStdForAllResponses(cyclenumber,:) = nanstd(Infos,1);
InfoSort = sort(Infos,1);
InfoLow(cyclenumber,:) = prctile(InfoSort,5,1);
InfoHigh(cyclenumber,:) = prctile(InfoSort,95,1);
end %samplecycle loop

% save the data from the information run
cd('/Users/frick/Documents/MATLAB')
save(strcat('INFORMATIONz',date,'-',num2str(dimensions),'-',num2str(iter),'d.mat'));
end



function DataCell = makePRofScellarray(DataCell,SigFeaturesArray,assembledFeatures)
%construct a cellarray of PRofS where dim1 is the signal number and dim2 is
%the feature of interest to analyze (such as cellular response)
for i = 1:size(DataCell,1)
    for j=1:length(assembledFeatures)
    DataI = SigFeaturesArray{i};
    CH = DataI(assembledFeatures{j},:);
    DataCell{i,j} = CH;
    end
end
end

function infoMatrix = calculateInfo(PRofScellarray,SampSize,numberOfSignals,PofSvector,lengthOfFeatures,k,DigK,entropyRgS,iter,assembledFeatures,dimensions)
infoMatrix = zeros(iter,lengthOfFeatures);




parfor i = 1:iter
% for i = 1:iter

uncx = struct();


    allIndices = cell(1,length(SampSize));
    sampleIndices = cell(1,length(SampSize));
    for cyc = 1:length(SampSize)
    allIndices{cyc} = 1:SampSize(cyc);
    sampleIndices{cyc} = false(1,SampSize(cyc));
    end

    
    Ni = min(cellfun(@length,allIndices,'UniformOutput',1));
    while Ni>0
        detSignal = datasample(1:numberOfSignals, 1, 'Weights', abs(PofSvector));
        vec = allIndices{detSignal}; %vec is a vector of indices for a given P(R|S=s)
        randIndex = datasample(vec, 1); %choose one index from this vector
        vec(vec==randIndex) = []; % then remove it (we are not resampling data)
        allIndices{detSignal} = vec; %save this new vector to allIndices

        sampvec = sampleIndices{detSignal}; %sampvec will become populated with ones at the location of the indices chosen from vec.
        sampvec(randIndex) = true; %set equal to 1
        sampleIndices{detSignal} = sampvec; %save this new vector to sampleIndices

        Ni = min(cellfun(@length,allIndices,'UniformOutput',1)); %determine the number of data indices remaining in each P(R|S=s)
    end
    

        
    for ooo = 1:lengthOfFeatures
        entropyrgsvector = entropyRgS(:,ooo);
        for d = 1:numberOfSignals
            PRofSgiven_s_vector=PRofScellarray{d,ooo};
            sampleIdx =  sampleIndices{d};
            PRofSgiven_s_vector(:,sampleIdx);
%             lengthofPRofSgsv = 1:size(PRofSgiven_s_vector,2);
%             sampleIndices = datasample(lengthofPRofSgsv,NewSamp(d),'Replace',false);
%             uncx(d).PRofSgiven_s_vector = PRofSgiven_s_vector(:,sampleIndices)';
            uncx(d).PRofSgiven_s_vector = PRofSgiven_s_vector(:,sampleIdx)';
        end
        
        Uncondx = vertcat(uncx.PRofSgiven_s_vector);
        UncondX = Uncondx';
        UncondResponsesGivenSvector =UncondX';
        
        if sum(sum(isnan(UncondResponsesGivenSvector)))>0
            idx = reshape(isnan(UncondX),size(UncondX));
            [~,col,~] = find(idx==1);
            [~,dim2] = size(UncondX);
            nline = 1:dim2;
            nline(col)=[];
            UncondResponsesGivenSvector = UncondX(:,nline)';  
        end
    
    
        if dimensions >1 %changed to 1 from 2 on 2015_12_07
            featuresSelected = assembledFeatures{ooo};
            dmnsn = determineDimension(featuresSelected,UncondResponsesGivenSvector);
        else
            dmnsn=1;
        end
        
        cd = (pi.^(dmnsn./2))./gamma(1+dmnsn./2)./2.^dmnsn; %volume of d-dimensional unit ball.

        if dimensions >1
            urgs = [];
            UncondResponsesGivenSvector = UncondResponsesGivenSvector;
        elseif assembledFeatures{ooo}==2 || assembledFeatures{ooo}==22;
            urgs = unique(UncondResponsesGivenSvector);
            UncondResponsesGivenSvector = urgs;
        else
            urgs = [];
            UncondResponsesGivenSvector = UncondResponsesGivenSvector;    
        end

        [~, Eps] = knnsearch(UncondResponsesGivenSvector,UncondResponsesGivenSvector,'K', k);
        Eps = Eps(:,k).*2; %Eps is the distance times two. 
        entropyR = (dmnsn./length(UncondResponsesGivenSvector))*sum(log(Eps)) - DigK + psi(length(UncondResponsesGivenSvector)) + log(cd);
        infoMatrix(i,ooo) = -(entropyR-sum(PofSvector'.*entropyrgsvector));  %sum(PofSvector'.*entropyrgsvector) is %step2[form the weighted average of h(R|S=s) using the the P(S), signal distribution]
        
        
        
    end %ooo loop
end  %iter loop
end



% % % SampSize = cellfun(@length,Data);
% % % NewSamp = zeros(size(SampSize));
% % % while all(SampSize-NewSamp)
% % % NewSamp = NewSamp + hist(datasample(1:numberOfSignals, min(SampSize-NewSamp), 'Weights', abs(SigProb)),1:numberOfSignals)';
% % % %take histogram of...(the following, 1:numberOfSignals)
% % % %The number of data is 1:N
% % % %The number of data values sampled is = min(SampSize-NewSamp)
% % % %'Weights' = Signal Probability
% % % %%%This should generate a new vector with numbers that will correspond to
% % % %%%how many data points need to be sampled corresponding to the signaling
% % % %%%probablity to generate the unconditional estimation. 
% % % end
% % % UncondX = [];
% % % for d = 1:numberOfSignals
% % %     UncondX = [UncondX datasample(Data{d}, NewSamp(d),'Replace', false)];
% % %     %sample a set number (given by NewSamp(d), which corresponds to the
% % %     %number.*weight for a signal probability from the data corresponding to
% % %     %that signal). 
% % % end
% % % figure, hist(UncondX,30)



function SamSigs = generateInputDistributions(scalarLength,smthness,number_of_PofS,lowsigma,highsigma)
Asigma = linspace(lowsigma,highsigma,number_of_PofS./2);
CDsigma = sort(Asigma,'descend');
Asigma = horzcat(Asigma,CDsigma); %lowsigma->highsigma->lowsigma with length of number_of_PofS
CDsigma = horzcat(CDsigma,ones(size(Asigma)).*lowsigma); %highsigma -> lowsigma, for half of length of number_of_PofS and then lowsigma->lowsigma for rest of lenght of number_of_PofS
%by altering the sigma we can get wide or narrow distributions so that when
%we run the script below, the dominant distributions will either be
%centered at 1 (unimodal), 0 and 2 (bimodal) or 0,1,and,2 (trimodal).
SamSigs = zeros(number_of_PofS,scalarLength);
    for i=1:length(Asigma)
        A = normrnd(1,Asigma(i),smthness,1);  %a vector with length#smoothness assembled by randomly sampling from a normaldistribution with mu=1 and sigma = Asigma(i)
        C = normrnd(0,CDsigma(i),smthness,1); %a vector with length#smoothness assembled by randomly sampling from a normaldistribution with mu=0 and sigma = CDsigma(i)
        D = normrnd(2,CDsigma(i),smthness,1); %a vector with length#smoothness assembled by randomly sampling from a normaldistribution with mu=2 and sigma = CDsigma(i)
        E = vertcat(A,C,D); %assembles A,C,D into a long vector
        h=histogram(E); 
        h.Normalization = 'probability';
        h.BinEdges = linspace(0,2,scalarLength+1);
        h.Normalization = 'probability';
        SamSigs(i,:) = h.Values; %take the values from the histogram 
    end
f=h.Parent;
ff=f.Parent;
ff.delete;
imagesc(SamSigs) %displays image of the assembled input distributions. These are P(S). where S is the signal (or condition)
end

function dmnsn = determineDimension(featuresSelected,responsesGivenSvector)
[~,epsdim] = knnsearch(featuresSelected',featuresSelected','K',2);    
    if epsdim == 0  
        dmnsn = 1;
    else
        epsdim = epsdim(:,2);
        dimd = sum(epsdim==0);
        if dimd ==3
            dmnsn = 1;
        elseif dimd ==2
            dmnsn = 2;
        else 
            dmnsn = size(responsesGivenSvector,2); %dimension is either 1D or 2D
        end
    end
end

function chosenScalars = chooseQFfunct(SCALAR,featuresVector,dimensions)
cycle = 1;
if isempty(featuresVector)   
    featuresVector = 1:size(SCALAR{1},1);   

    if dimensions ==1
        chosenScalars = cell(1,length(featuresVector).^dimensions);
        for i=featuresVector
        chosenScalars{cycle} = [i];
        cycle=cycle+1;
        end  
    end
    
    if dimensions ==2
        chosenScalars = cell(1,length(featuresVector).^dimensions);
        for i=featuresVector
            for j =featuresVector
            chosenScalars{cycle} = [i j];
            cycle=cycle+1;
            end
        end
    end

    if dimensions ==3
    chosenScalars = cell(1,length(featuresVector).^dimensions);
        for i=featuresVector
            for j =featuresVector
                for kk =featuresVector
                chosenScalars{cycle} = [i j kk];
                cycle=cycle+1;
                end
            end
        end 
    end
      

else

    if dimensions ==1
        chosenScalars = cell(1,length(featuresVector).^dimensions);
        for i=featuresVector
        chosenScalars{cycle} = [i];
        cycle=cycle+1;
        end
    end


    if dimensions ==2
        chosenScalars = cell(1,length(featuresVector).^dimensions);
        for i=featuresVector
            for j =featuresVector
            chosenScalars{cycle} = [i j];
            cycle=cycle+1;
            end
        end
    end

    if dimensions ==3
        chosenScalars = cell(1,length(featuresVector).^dimensions);
        for i=featuresVector
            for j =featuresVector
                for kk =featuresVector
                chosenScalars{cycle} = [i j kk];
                cycle=cycle+1;
                end
            end
        end
    end


end
end