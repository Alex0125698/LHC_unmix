% Hyperspectral Generator & Unmixing
% by A.S. Woodcock; 22/AUG/19

clear variables
close all
clc

%%% Load Data

Y = csvread('stopScaled_alls.csv');

hist(Y(:,8),30);
% hold on
% hist(Y(:,2),60);
% hist(Y(:,3),60);
% hold off
% Y = Y';
% save('stopScaled_alls.mat','Y');

% load('Indian_pines.mat');
% d = indian_pines / max(indian_pines,[],'all');
% Y2 = reshape(d,size(d,1)*size(d,2),size(d,3));

%%% Pixel Generation

Npixels  = 1000;
Nsigs = 4;
Nfreqs = 22;

[x,u] = generateSpectra(Nsigs,Nfreqs);

% Mixing matrix Nfreqs x Nsigs
M = u';

[Y,S,W] = PixelGenerator(M,Npixels,1,1);

% Write ENVI files
% data=reshape(Y',10,Npixels/10,Nfreqs);
% info=enviinfo(data);
% info.wavelength_units = 'None';
% info.z_plot_titles = 'None, Magnitude';
% dummyNames = sprintf('%.0f, ', (0:Nfreqs-1));
% dummyNames = strcat('{',dummyNames(1:end-2),'}');
% info.band_names = dummyNames;
% info.wavelength = dummyNames;
% enviwrite(data,info,'test.dat');

% Write matlab data
save('easy.mat','Y','M','S');

[SNR_SD,xi] = PCA(Y,W);

figure
semilogy(xi,SNR_SD,'-',xi,SNR_SD,'O');

it = logspace(2,4,10);
error = zeros(10,1);

for i = 1:10
    [M2,S2] = LS_NMF(Y,Nsigs,it(i));
    error(i) = sum((M*S - M2*S2).^2,'all');
end

% disp(error);

M3 = N_FINDER(Y);

figure
subplot(1,2,1);
plot(x',u','-',x',u','O');
title('Spectral Signatures');
xlabel('rel. frequency');
ylabel('rel. strength');
subplot(1,2,2);
scatter3(Y(1,:),Y(2,:),Y(3,:),'.');
hold on
scatter3(M(1,:),M(2,:),M(3,:),'x');
% hold on
% scatter3(M2(1,:),M2(2,:),M2(3,:),36);
hold on
scatter3(M3(1,:),M3(2,:),M3(3,:),36);
hold off
title('Pixel Data (First 3 freq. Components)')
xlabel('Freq 1');
ylabel('Freq 2');
zlabel('Freq 3');
legend('pixels','endmembers','LS-NMF','N-FINDR');

%%% Principle Component Analysis

function [SNR_SD,xi] = PCA(Y,W)
    % nFreqs x nPixels -> n2 x nPixels
    A = Y';
%     A = A - mean(A);
%     A = A ./ max(abs(A));
    Ry = cov(A);
    Rw = cov(W');
    [V,D] = eig(Ry);
    d = diag(D);
    
    nVecs = size(Ry,1);
    SNR_SD = zeros(1,nVecs);
    
    for i=1:nVecs        
        SNR_SD(i) = d(i) / (V(:,i)' * Rw * V(:,i));        
    end
    
    xi = 1:nVecs;
      
end

%%% Constrained Lee and Seung's multiplicative update rule

function [W,H] = LS_NMF(V,innerDim, iterations)
    % V = W x H
    
    n = size(V,1);
    m = size(V,2);
    
    W = rand(n, innerDim);
    H = ones(innerDim, m);
    H = H / sum(H);
    % H nxi
    % W ixm
    % V nxm
    
    for j=1:iterations
        H = H .* (W' * V) ./ ((W' * W) * H);
        W = W .* (V * H') ./ (W * (H * H'));
        H = H*0.99 + (H / sum(H))*0.01;
    end

end

%%% N-FINDR Algorithm

function M = N_FINDER(Y)
    % must have nFreqs = nSigs - 1
    
    p = size(Y,2);
    n = size(Y,1);
    e = 1:n+1;
    
    Vold = 0;
    for k=1:3
        for i=1:n+1
            for j=1:p

                e2 = e;
                e2(i) = j;

                V = SimplexVolume(Y(:,e2));

                if V > Vold
                    Vold = V;
                    e = e2;                
                end
            end
        end
    end
    
    M = Y(:,e);
end

function V = SimplexVolume(P)
    % P is a matrix whose columns are the points of an
    % n-dimensional simplex (there must be n+1 points)
    % Each row represents a dimension and they can be 
    % placed in any order

    n = size(P,2);
    V = abs(det([ones(1,n);P])) / factorial(n-1);

end
