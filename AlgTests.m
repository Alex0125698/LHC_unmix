% Runs a vairety of unmixing Algorithms on a variety of
% Hyperspectral Images and plots results
%
% By A.S. Woodcock (03/NOV/2019)

clearvars
close all
clc

addpath('data');
addpath('algs');

% number of hyperspectral images
nImages = 1;
% number of unmixing algorithms
nAlgs = 7;
% performances
perf = zeros(nAlgs,nImages);
% name of each unmixing algorithm
algNames = {'DUSAL','DECA','VCA','NFINDR','MVSA','SISAL','BSDMM'};
% name of each hyperspectral image
imgNames = {'Easy','Standard','Noisy','Mixed','Lopsided','Hard'};
% imgNames = {'LHC-1D Hist'};



% load('imgDemoNoise');
% perf(:,1) = runAlgs(zeros(3,4),Y);
% load('imgHSDemo');
% perf(:,2) = runAlgs(zeros(3,4),Y);

load('imgEasy');
% Ek = ne_hysime(Y)
perf(:,1) = runAlgs(M,Y);
load('imgStandard');
% Ek = ne_hysime(Y)
perf(:,2) = runAlgs(M,Y);
load('imgNoisy');
% Ek = ne_hysime(Y)
perf(:,3) = runAlgs(M,Y);
load('imgMixed');
% Ek = ne_hysime(Y)
perf(:,4) = runAlgs(M,Y);
load('imgLopsided');
% Ek = ne_hysime(Y)
perf(:,5) = runAlgs(M,Y);
load('imgHard');
% Ek = ne_hysime(Y)
perf(:,6) = runAlgs(M,Y);

% load('EventData1DHist');
% Ek = ne_hysime(Y)

% perf(:,1) = runAlgs(M,Y);

figure
% c = categorical(algNames);
bar(categorical(algNames),perf);
legend(imgNames);
% ylabel('Performance');
title('Unmixing Performance (0=bad,100=perfect)');
set(gca, 'YGrid', 'on', 'XGrid', 'off')
set(gca,'ytick',linspace(0,100,11))
ylim([0,100]);
% grid on

function perf = runAlgs(M,Y)

    [nBands,nEndmems] = size(M);
    nAlgs = 7;
    M_algs = zeros(nAlgs,nBands,nEndmems);

    Deca_Parameters = struct('Endmembers',nEndmems, ... 
                            'Tol_th',1e-5, ...
                            'Number_iteration_max',40, ...
                            'Display_figure','off', ...
                            'Verbose','off'...
                            );

    Malg = hu_dusal(Y,nEndmems);
    M_algs(1,:,:) = reorder(M,Malg);
%     Malg = hu_deca_mdl(Y,Deca_Parameters);
    Malg = hu_dusal(Y,nEndmems);
    M_algs(2,:,:) = reorder(M,Malg);
    Malg = hu_vca(Y,'Endmembers',nEndmems);
    M_algs(3,:,:) = reorder(M,Malg);
    Malg = hu_nfindr(Y,nEndmems);
    M_algs(4,:,:) = reorder(M,Malg);
    Malg = mvsa(Y,nEndmems);
    M_algs(5,:,:) = reorder(M,Malg);
    Malg = hu_sisal(Y,nEndmems);
    M_algs(6,:,:) = reorder(M,Malg);
    % Malg = nmf_als(Y,nEndmems,100,0);
    % M_algs(7,:,:) = reorder(M,Malg);
    
    % our final algorithm is in Python
    % save data that our python script will load
    DataGen.saveHSImage('hsuDataIn',0,Y,M,0,0);
    % run the script and wait until it completes
    system('python unmixing.py');
    % load the results
    load('data/hsuDataOut','Malg');
    M_algs(7,:,:) = reorder(M,Malg);
        
%     [~,I] = sort(sum(M_algs,2),3);
%     for i=1:nAlgs
%         M_algs(i,:,:) = M_algs(i,:,squeeze(I(i,1,:)));
%     end
%     [~,I] = sort(sum(M));
%     M = M(:,I);

    figure
    subplot(1,3,1);
    plotEndmembers(M,squeeze(M_algs(1,:,:)),'DUSAL')

    subplot(1,3,2);
    plotEndmembers(M,squeeze(M_algs(3,:,:)),'NFINDR')
    
    guess = mean(Y,2) + zeros(size(M));
    
    subplot(1,3,3);
    plotEndmembers(M,guess,'Guess')

    err_guess = calcPerformance(M,guess);
    
    perf = zeros(nAlgs,1);
    
    for i=1:nAlgs
        err = calcPerformance(M,M_algs(i,:,:));
        perf(i) = 100 * (1 - err/err_guess);
        if perf(i) < 0, perf(i) = 0; end
    end

end

function plotEndmembers(M,M_alg,algName)
    nBands = size(M,1);
    hold on
    plot(1:nBands,M_alg);
    ax = gca;
    ax.ColorOrderIndex = 1;
    plot(1:nBands,M,'--');
    legend([plot(nan,nan,'-k'),plot(nan,nan,'--k')], algName, 'True');
    title(strcat(algName,' (Endmember Signatures)'));
    xlabel('Frequency');
    ylabel('Magnitude');
    hold off
end

function err = calcPerformance(M,M_alg)
    M = M(:);
    M_alg = M_alg(:);
    err = sum(abs(M-M_alg)) / length(M);
end

function M_alg = reorder(M,M_alg)
    
    [nBands,nEndmems] = size(M);
    
    orders = zeros(nEndmems,1);
    
    for i=1:nEndmems
        e_alg = M_alg(:,i);
        err = Inf;
        orderI = nan;
        
        for j=1:nEndmems
            if sum(orders == j) ~= 0, continue, end
            
            ej = M(:,j);
            err2 = sum(abs(ej-e_alg));
            
            if err2 < err
                orderI=j; 
                err = err2; 
            end         
        end
        orders(i) = orderI;
    end
    
%     M_alg = M_alg(:,orders);
    M_alg(:,orders) = M_alg;
end