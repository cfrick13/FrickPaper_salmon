function cmap = cmapGenerator(lengthofdata)

cno=[];
if isEven(lengthofdata) == 0
    lengthofdata = lengthofdata+1;
end

dim = [lengthofdata 3];
lumiRange = [85 85 1]; % luminance
alphaRange = [50 -50 0]; %green to red axis
betaRange = [50 50 0]; % blue to yellow axis
labalt = zeros(dim);
labalt(:,1) = linspace(lumiRange(1),lumiRange(2),dim(1));
labalt(:,2) = linspace(alphaRange(1),alphaRange(2),dim(1));
labalt(:,3) = linspace(betaRange(1),betaRange(2),dim(1));

if lumiRange(3) == 1
    labalt(:,1) = flipud(labalt(:,1));
end
if alphaRange(3) == 1
    labalt(:,2) = flipud(labalt(:,2));
end
if betaRange(3) == 1
    labalt(:,3) = flipud(labalt(:,3));
end

% cmapfinaltwo = lab2rgb(labalt,'ColorSpace','adobe-rgb-1998');
cmapfinaltwo = lab2rgb(labalt);
% cmapimage = reshape(cmapfinaltwo,[1 size(cmapfinaltwo,1) size(cmapfinaltwo,2)]);
% figure
% imshow(vertcat(cmapimage,cmapimage,cmapimage,cmapimage,cmapimage,cmapimage,cmapimage,cmapimage));
cmapfinal = cmapfinaltwo;
cmapfinal(cmapfinal>1)=1;
cmapfinal(cmapfinal<0)=0;
cmap = cmapfinal;

cmap = flipud(cmap);
%% gray - blue - green

if lengthofdata ==4
   lengthofdata = 6; 
   cno=1;
end

dim = [lengthofdata 3];
lumiRange = [50 70 0]; % luminance
alphaRange = [0 5 -80 0]; %green to red axis
betaRange = [0 -50 20 0]; % blue to yellow axis
labalt = zeros(dim);
labalt(:,1) = linspace(lumiRange(1),lumiRange(2),dim(1));
labalt(:,2) = horzcat(linspace(alphaRange(1),alphaRange(2),dim(1)./2),linspace(alphaRange(2),alphaRange(3),dim(1)./2));
% labalt(:,3) = linspace(betaRange(1),betaRange(2),dim(1));
labalt(:,3) = horzcat(linspace(betaRange(1),betaRange(2),dim(1)./2),linspace(betaRange(2),betaRange(3),dim(1)./2));

if lumiRange(3) == 1
    labalt(:,1) = flipud(labalt(:,1));
end
if alphaRange(3) == 1
    labalt(:,2) = flipud(labalt(:,2));
end
if betaRange(3) == 1
    labalt(:,3) = flipud(labalt(:,3));
end

% cmapfinaltwo = lab2rgb(labalt,'ColorSpace','adobe-rgb-1998');
cmapfinaltwo = lab2rgb(labalt);
% cmapimage = reshape(cmapfinaltwo,[1 size(cmapfinaltwo,1) size(cmapfinaltwo,2)]);
% figure
% imshow(vertcat(cmapimage,cmapimage,cmapimage,cmapimage,cmapimage,cmapimage,cmapimage,cmapimage));
cmapfinal = cmapfinaltwo;
cmapfinal(cmapfinal>1)=1;
cmapfinal(cmapfinal<0)=0;
cmap = cmapfinal;

cmap = flipud(cmap);

if ~isempty(cno)
    cmap = cmap([1 2 4 6],:);
end

%% gray - orange - green
% dim = [lengthofdata 3];
% lumiRange = [65 65 1]; % luminance
% alphaRange = [0 20 -60 0]; %green to red axis
% betaRange = [0 80 10 0]; % blue to yellow axis
% labalt = zeros(dim);
% labalt(:,1) = linspace(lumiRange(1),lumiRange(2),dim(1));
% labalt(:,2) = horzcat(linspace(alphaRange(1),alphaRange(2),dim(1)./2),linspace(alphaRange(2),alphaRange(3),dim(1)./2));
% % labalt(:,3) = linspace(betaRange(1),betaRange(2),dim(1));
% labalt(:,3) = horzcat(linspace(betaRange(1),betaRange(2),dim(1)./2),linspace(betaRange(2),betaRange(3),dim(1)./2));
% 
% if lumiRange(3) == 1
%     labalt(:,1) = flipud(labalt(:,1));
% end
% if alphaRange(3) == 1
%     labalt(:,2) = flipud(labalt(:,2));
% end
% if betaRange(3) == 1
%     labalt(:,3) = flipud(labalt(:,3));
% end
% 
% % cmapfinaltwo = lab2rgb(labalt,'ColorSpace','adobe-rgb-1998');
% cmapfinaltwo = lab2rgb(labalt);
% % cmapimage = reshape(cmapfinaltwo,[1 size(cmapfinaltwo,1) size(cmapfinaltwo,2)]);
% % figure
% % imshow(vertcat(cmapimage,cmapimage,cmapimage,cmapimage,cmapimage,cmapimage,cmapimage,cmapimage));
% cmapfinal = cmapfinaltwo;
% cmapfinal(cmapfinal>1)=1;
% cmapfinal(cmapfinal<0)=0;
% cmap = cmapfinal;
% 
% cmap = flipud(cmap);
end
