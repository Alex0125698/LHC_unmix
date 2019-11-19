% Generates a variety of test data for Hyperspectral Unmixing Algorithms
%
% Author: A.S. Woodcock 22/OCT/2019

% dummy class to allow export of multiple functions
% run using: DataGen.functionName(...)
classdef DataGen
methods(Static)
    
    
function run
    % Generates + saves data
    
    clearvars
    close all
    clc

    nBands = 20; % number of frequency bands
    nEndmems = 4; % number of endmembers/signatures
    height = 50; % height of image
    width = 50; % width of image

    M = DataGen.genEndmems(nBands,nEndmems);

    figure
    hold on
    plot(1:nBands,M,'-');
    ax = gca;
    ax.ColorOrderIndex = 1;
    plot(1:nBands,M,'o');
    title('Endmember Frequency signatures');
    xlabel('Frequency');
    ylabel('Magnitude');
    hold off

    savefig('data/imgAllSigs');
    saveas(gcf,'data/imgAllSigs.png');



    % Generate and save a variety of hyperspectral Images

    figure('Name', 'Generated HSI Data');

    [HSI,Y,S,W] = DataGen.genHSImage(M,width,height,3,0,0);
    DataGen.saveHSImage('imgEasy',HSI,Y,M,S,W);
    subplot(2,3,1);
    DataGen.plotData(Y,M,'Easy HSI');

    [HSI,Y,S,W] = DataGen.genHSImage(M,width,height,+0.5,0.5,0.4);
    DataGen.saveHSImage('imgStandard',HSI,Y,M,S,W);
    subplot(2,3,2);
    DataGen.plotData(Y,M,'Standard HSI');

    [HSI,Y,S,W] = DataGen.genHSImage(M,width,height,+0.5,6,0);
    DataGen.saveHSImage('imgNoisy',HSI,Y,M,S,W);
    subplot(2,3,3);
    DataGen.plotData(Y,M,'Noisy HSI');

    [HSI,Y,S,W] = DataGen.genHSImage(M,width,height,-11,1,0.2);
    DataGen.saveHSImage('imgMixed',HSI,Y,M,S,W);
    subplot(2,3,4);
    DataGen.plotData(Y,M,'Well-mixed HSI');

    [HSI,Y,S,W] = DataGen.genHSImage(M,width,height,0,0.5,6);
    DataGen.saveHSImage('imgLopsided',HSI,Y,M,S,W);
    subplot(2,3,5);
    DataGen.plotData(Y,M,'Lopsided HSI');

    [HSI,Y,S,W] = DataGen.genHSImage(M,width,height,-11,5.5,12);
    DataGen.saveHSImage('imgHard',HSI,Y,M,S,W);
    subplot(2,3,6);
    DataGen.plotData(Y,M,'Hard HSI');

    img = imread('imgHSdemo.png');
    Y = reshape(img,size(img,1)*size(img,2),size(img,3))';
    DataGen.saveHSImage('imgHSdemo',HSI,Y,0,0,0);
    % subplot(2,3,5);
    % plotData(Y,[],'Img Demo with endmembers');

    img = imread('imgDemoNoise.png');
    Y = reshape(img,size(img,1)*size(img,2),size(img,3))';
    DataGen.saveHSImage('imgDemoNoise',HSI,Y,0,0,0);
    % subplot(2,3,6);
    % DataGen.plotData(Y,[],'Img Demo random');

    savefig('data/imgAllData');
    saveas(gcf,'data/imgAllData.png');

end

function imageExamples
    
    nBands = 3; % number of frequency bands
    nEndmems = 4; % number of endmembers/signatures
    height = 30; % height of image
    width = 30; % width of image

    M = DataGen.genEndmems(nBands,nEndmems);
    
	[HSI,Y,S,W] = DataGen.genHSImage(M,width,height,3,0,0); 
    figure
    subplot(1,3,1)
    DataGen.plotData(Y,M,'Linear HSI');
    
    subplot(1,3,2)
    Y = rand(3,10000)*255;
    scatter3(Y(1,:),Y(2,:),Y(3,:),'.');
    title('Normal Image');
    axis([-50 300 -50 300 0 255]);
    xlabel('band 1');
    ylabel('band 2');
    zlabel('band 3');
    
    subplot(1,3,3)
    hold on
    Y1 = (randn(2,2000)*0.9 + 0.1);
    Y2 = (randn(2,1000)*0.7 + 0.8);
    Y3 = (randn(2,300)*0.3 - 0.4);
    scatter(Y1(1,:),Y1(2,:),'.');
    scatter(Y2(1,:),Y2(2,:),'.');
    scatter(Y3(1,:),Y3(2,:),'.');
    title('Image of gaussian distributions');
%     axis([-50 300 -50 300 0 255]);
    xlabel('band 1');
    ylabel('band 2');
    hold off
end

function saveColliderData()
    names = ["lep1pt", "jet1pt", "jet2pt", "met", "jet1phi",...
         "jet2phi", "mtbmin", "mtbmax", "drbb", "mt2a",...
         "mlbmin", "mlbmax", "ht", "metphi", "leptonIsElectro n"];
     
    data = csvread('stopScaled_alls.csv');
    save('data/stop_signal','names','data');
    data = csvread('topScaled_allb.csv');
    save('data/top_background','names','data');
end

function plotData(Y,M,name)

    hold on
    scatter3(Y(1,:),Y(2,:),Y(3,:), 'b.');
    if ~isempty(M)
        x = M(1,:);
        y = M(2,:);
        z = M(3,:);
        scatter3(M(1,:),M(2,:),M(3,:),'rX');
        indices = [1;2;3;4;1;3;2;4];
        plot3(x(indices),y(indices),z(indices),'k-');
    end
    title(name);
    xlabel('band 1');
    ylabel('band 2');
    zlabel('band 3');
    hold off

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

function [HSI,Y,S,W] = genHSImage(M,width,height,purity,noise,lopsided)
    % Generate an random hyperspectral image using
    % a known set of engmember signatures
    % ---- Input ----
    % M is the 'Mixing matrix' (nBands x nEndmems)
    % height is the length of the first dim of the image
    % width is the length of the second dim of the image
    % purity controls how pure the pixels are:
    %   less-than 0 -> more mixed,
    %   greater-than 0 -> pure pixels,
    %   equal to 0 -> 
    %   recommended range: -5 to +5
    % Noise controls how noisy the pixels are (gaussian)
    %   0 -> no noise; g.t. 0 -> noisy; range 0 to 10
    % lopsided controls how the data is distributed between endmembers
    %   0 -> no endmembers favoured
    %   greater than 0 -> data is often closer to particular endmembers
    %   recommended range: 0 to 10
    % ---- Output ----
    % Y is the final image matrix (nBands x nPixels)
    % S is the abundance matrix (nEndmems x nPixels)
    % W is the noise matrix (nBands x nPixels)
    %
    % Author: A.S. Woodcock 22/OCT/2019
    
    [nBands,nEndmems] = size(M);
    nPixels = width * height;
    
    % we divide the interval [0,1] into nEndmems subintervals
    % at random;
    dividers = rand(nEndmems-1,nPixels);
    
    if purity >= 0
        % make pixels more pure rather than random
        dividers = dividers.^(1+purity);
    end
    
    dividers = sort(dividers);
    x = [dividers;ones(1,nPixels)];
    y = [zeros(1,nPixels);dividers];
    S = x - y;
    
    % permutate the data so that there is no bias towards
    % the last endmember
    for i = 1:size(S,2)       
        S(:,i) = S(randperm(nEndmems),i);       
    end
    
    % add bias towards endmembers
    S = S .* (exp(-(lopsided+2) * 0.3 * (0:nEndmems-1)')*4*lopsided+1);
    S = S ./ sum(S);
    
    % make pixel data less pure
    if purity < 0
        S = S - 0.3*purity*rand(nEndmems,nPixels) + abs(purity)/20;
        S = S ./ sum(S);
    end
    S = S ./ sum(S);

    % add gaussian noise
    W = randn(nBands,nPixels) * 0.005 * noise;
    Y = M*S + W;
    HSI = reshape(Y',width,height,nBands);
end

function M = genEndmems(nBands,nEndmems)
    % Generate an random 'Mixing matrix' containting 
    % endmember signatures; size = [nBands,nEndmems]
    % ---- Input ----
    % nBands = number of frequency components/bands
    % nEndmems = number of endmembers/signatures
    % ---- Output ----
    % M = the (nBands x nEndmems) Mixing matrix
    % 
    % Author: A.S. Woodcock 22/OCT/2019

    M = zeros(nBands,nEndmems);
    
    for i=1:nEndmems
        % some complicated parameters to give us random looking endmembers
        nComps = rand(1)*15 + 22;  
        freqs = i + 0.6*(0.05*i/nComps + 1.13).^(1:nComps);
        mags = (flip(1:nComps)/nComps).^(0.6+0.05*i/nComps);
        mags = mags * (rand(1)+0.01)^3;
        
        M(:,i) = DataGen.noiseOctave(nBands,freqs,mags)';
    end
    
    M = (M - min(M));
    M = M ./ max(M);

    M = M + rand(1,nEndmems)*0.4;
    M = M .* (rand(1,nEndmems)+0.5);
    
end

function u = noiseOctave(nPoints,freqs,mags)
    % Generate a vector of smooth'ish noise 
    % with a number of rough freq. components each with 
    % a set magnitude
    % inspired by Perlin noise
    % ---- Input ----
    % nPoints = number of points to generate
    % freq = rough digital frequency of noise
    % ---- Output ----
    % u = the vector of smooth noise
    % 
    % Author: A.S. Woodcock 22/OCT/2019
    
    u = zeros(nPoints,1);
    
    for i=1:length(freqs)
        u = u + DataGen.smoothNoise1D(nPoints,freqs(i)) * mags(i);
    end
end

function u = smoothNoise1D(nPoints,freq)
    % Generate a vector of smooth'ish noise which
    % varies roughly at a set frequency
    % inspired by Perlin noise
    % ---- Input ----
    % nPoints = number of points to generate
    % freq = rough digital frequency of noise
    % ---- Output ----
    % u = the vector of smooth noise (-1 to +1)
    % 
    % Author: A.S. Woodcock 22/OCT/2019

    fade = @(x) x.*x.*x.*(x.*(x.*6-15)+10);
    
    r = rand(round(freq)+2,1);
    x = (0:nPoints-1)';
    
    index = floor(x*freq/(nPoints-1));
    grad1 = r(index+1);
    grad2 = r(index+2);
        
    pos = mod(x*freq/(nPoints-1),1);
    v1 = fade(pos);
    v2 = fade(1-pos);
    
    u = grad1.*v2 + grad2.*v1;
end


end
end