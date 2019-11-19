clearvars
close all
clc

addpath('data');
addpath('algs');

load('stop_signal');
dataS = data(:,1:14);
load('top_background');
dataB = data(:,1:14);

nVars = size(dataB,2);

names = ["p_T^l^e^p^1","p_T^j^1","p_T^j^2","E_T^m^i^s^s","\phi^j^1","\phi^j^2","m_T^b^,^m^i^n",...
    "m_T^b^,^m^a^x","\DeltaR","m_T_2","m_b_l^m^i^n","m_b_l^m^a^x","H_T","\phi^E^m^i^s^s",...
    "leptonIsElectro n"];

% make sure events are in the same range for signal & background
% otherwise generated histograms are different
vMin = max([min(dataS);min(dataB)]);
vMax = min ([max(dataS);max(dataB)]);

% dataS = dataS(sum(dataS >= vMin,2) == nVars,:);
% dataS = dataS(sum(dataS <= vMax,2) == nVars,:);
% dataB = dataB(sum(dataB >= vMin,2) == nVars,:);
% dataB = dataB(sum(dataB <= vMax,2) == nVars,:);

[nSig, nVars] = size(dataS);
nBack = size(dataB,1);

nS = length(dataS);
nB = length(dataB);

nMin = min(nS,nB);

dataS = dataS(1:nMin,:);
dataB = dataB(1:nMin,:);

data = [dataS;dataB];

M = [mean(dataS);mean(dataB)]';
Y = [dataS;dataB]';

drawHistograms1D(dataS,dataB,names, vMin, vMax);
for i=1:nVars
    drawHistograms2D(dataS,dataB,names,i)
end

% createHistograms1D(dataS,dataB,names, vMin, vMax);
% 
% saveHSImage('LHCplain',0,Y,M,0,0);

function createHistograms1D(dataS,dataB,names, vMin, vMax)
    % represents LHC event variables as an image
    % using bins taken from 1D histograms of each
    % variable as the frequency bands

    % the variables which show differing behaviour between
    % signal and background
    interestingVars = [1 2 3 4 7 8 10 11 12];
    
    dataS = dataS(:,interestingVars);
    dataB = dataB(:,interestingVars);
    
    data = [dataS;dataB]; % combined signal + background data
    [nEvents,nVars] = size(data);
    
    % shuffle data
    data = data(randperm(nEvents),:);
    
    nHistPoints = 400; % number of events to use per histogram
    nHists = 40000; % number of histograms (of random points)
    nBins = 16; % number of histogram bins
    
    nBands = nBins*nVars; % number of image frequency bands
    nPixels = nHists; % number of image pixels
    
    Y = zeros(nBands,nPixels); % pixel matrix
    
    for i = 1:nHists
        tmp = randperm(nEvents);
        points = tmp(1:nHistPoints);
        
        for j = 1:nVars
            bandRange = ((j-1)*nBins+1):(j*nBins);
            Y(bandRange,i) = hist2(data(points,j), nBins, vMin(j), vMax(j))' / nHistPoints;
        end       
    end
    
    % Calculate true endmembers  
    nEndmems = 2; % number of endmembers
    M = zeros(nBands,nEndmems); % mixing matrix
        
    for j = 1:nVars
        bandRange = ((j-1)*nBins+1):(j*nBins);
        M(bandRange,1) = hist2(dataB(:,j),nBins, vMin(j), vMax(j))' / size(dataB,1);
        M(bandRange,2) = hist2(dataS(:,j),nBins, vMin(j), vMax(j))' / size(dataS,1);
    end
    
    saveHSImage('EventData1DHist',0,Y,M,0,0);
end

function drawHistograms2D(dataS,dataB,names,var1)
 
    figure
    nBins = 12;
    nVars = size(dataB,2);
    nFigRows = 3;
    nFigCols = ceil((nVars-1)/nFigRows);
%     levels = [0.01 0.03 0.1 0.3 0.55 0.9];
    levels = logspace(-2,0,9);
    
    for i=1:nVars
       var2 = i;
       if var1 == var2, continue, end
           
       subplot(nFigRows,nFigCols,i);
       hold on;
       N1 = hist3(dataB(:,[var1,var2]),'Nbins',[nBins,nBins]);
       N2 = hist3(dataS(:,[var1,var2]),'Nbins',[nBins,nBins]);
       N1 = N1 / size(dataB,1);
       N2 = N2 / size(dataS,1);
       
       m = max(max(max(N1)),max(max(N2)));
       N1 = N1 / m;
       N2 = N2 / m;
%        [i1,i2] = clip(counts1 / max(counts1),0.01);
%        [i3,i4] = clip(counts2 / max(counts2),0.01);       
%        counts1 = counts1(i1:i2);
%        counts2 = counts2(i3:i4);
%        centres1 = centres1(i1:i2);
%        centres2 = centres2(i3:i4);
       contour(N1,levels,'c');
       contour(N2,levels,'b');
       ylabel(names(var1));
       xlabel(names(var2));
%        title(names(i));
       hold off;
    end
end

function drawHistograms1D(dataS,dataB,names, vMin, vMax)

    figure
    nBins = 20;
    nVars = size(dataB,2);
    nFigRows = 4;
    nFigCols = ceil(nVars/nFigRows);
    
    for i=1:nVars
       subplot(nFigRows,nFigCols,i);
       hold on;
       [counts1,centres1] = hist2(dataB(:,i),nBins, vMin(i), vMax(i));
       [counts2,centres2] = hist2(dataS(:,i),nBins, vMin(i), vMax(i));       
       counts1 = counts1 / size(dataB,1);
       counts2 = counts2 / size(dataS,1);      
       [i1,i2] = clip(counts1 / max(counts1),-1);
       [i3,i4] = clip(counts2 / max(counts2),-1);       
       counts1 = counts1(i1:i2);
       counts2 = counts2(i3:i4);
       centres1 = centres1(i1:i2);
       centres2 = centres2(i3:i4);
       plot(centres1,counts1);
       plot(centres2,counts2); 
       title(names(i));
       hold off;
    end

end

function saveHSImage(name,HSI,Y,M,S,W)

    explanation = "Hyperspectral Image data:\n" + ...
        "\tHSI = Full hysperspectral image (width x height x nBands)\n" + ...
        "\tY = 'Flat' Hyperspectral image (nBands x nPixels)\n" + ...
        "\tM = Mixing Matrix (nBands x nEndmems)\n" + ...
        "\tS = Abundance Matrix (nEndmems x nPixels)\n" + ...
        "\tW = Noise Matrix (nBands x nPixels)\n" + ...
        "Note that: Y = MS + W\n" + ...
        "If data isn't avaliable, it is set to 0";
    
    [width,height,nBands] = size(HSI);
    [nEndmems,nPixels] = size(S);
    
    save(strcat('data/',name),'HSI','Y','M','S','W','width','height','nBands','nEndmems' ...
        ,'nEndmems','nPixels', 'explanation');

end

function bins = smoothHist(x,grid)
    
    grid = grid(:);
    x = x(:);
    n = length(x);
    stdev = 10/sqrt(n);
    bins = sum(exp(-((grid'-x)./stdev).^2))';

end

function [bins,centres] = hist2(x, nBins, xMin, xMax)

    x = x(:);    
    spacing = (xMax - xMin) / nBins;    
    binStart = spacing * (0:nBins-1) + xMin;
    binEnd = spacing * (1:nBins) + xMin;   
    bins = sum((x >= binStart) & (x < binEnd));
    centres = (binStart + binEnd)/2;
end

function [i1,i2] = clip(u,cutoff)

    i1 = length(u)*0+1;
    i2 = length(u);
    
    for i = i1:i2
        if u(i) > cutoff
            i1 = i;
            break;
        end
    end
    
    for i = flip(i1:i2)
        if u(i) > cutoff
            i2 = i;
            break;
        end
    end
end
