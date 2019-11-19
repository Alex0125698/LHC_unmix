clearvars
close all
clc

% classdef nDimEstimate
%     
% methods(Static)
% 
% function [x] = a
% x = 1;
% end
% 
% function [x] = b
% x = 2;
% end
% 
% function [x] = c(var)
% x = 3*var;
% end
% 
% end
% end

load('imgEasy');

% nDims = HFC(HSI,0.0001)
% nDims = HFC(HSI,0.001)
% nDims = HFC(HSI,0.01)
% nDims = HFC(HSI,0.05)
% nDims = HFC(HSI,0.1)
% nDims = HFC(HSI,0.5)
% nDims = HFC(HSI,1)

[endmemberindex,duration] = PPI(HSI,nEndmems);

M_est = nan(size(M));

for i=1:nEndmems
    M_est(:,i) = HSI(endmemberindex(i,1),endmemberindex(i,2),:);
end

figure
subplot(1,2,1);
plot(1:nBands,M);
title('True endmember signatures');
subplot(1,2,2);
plot(1:nBands,M_est);
title('Estimated endmember signatures');


% Purpose: Estimate number of endmembers
% Algorithm name: Harsanyi–Farrand–Chang (HFC) method
% Sensitivity to noise: moderate
function nDims = HFC(imgXYN,t)
    %
    % HFC gives the VD number estimated by given false alarm property using HFC
    % method.
    %
    % There are two parameters,HFC(imgXYN,t) where the imgXYN is the
    % Hyperspectral image cube, which is a 3-D data matrix
    % [XX,YY,bnd] = size(imgXYN), XX YY are the image size,
    % bnd is the band number of the image cube.
    % t is the false alarm probability.
    %
    % HFC uses the HFC algorithm developed by Dr. Chein-I Chang,
    % see http://www.umbc.edu/rssipl/. The Matlab code was
    % programmed by Jing Wang in Remote Sensing Signal and
    % Image Processing Lab.
    %
    [nX,nY,nBands] = size(imgXYN);
    nPixels = nX*nY;
    r = (reshape(imgXYN,nPixels,nBands))';
    R = (r*r')/nPixels;
    u = mean(r,2);
    K = R-u*u';
    
    %======HFC=====
    D1=sort(eig(R));
    D2=sort(eig(K));
    sita=sqrt((D1.^2+D2.^2)*2/nPixels);
    P_fa=t;
    Threshold=(sqrt(2))*sita*erfinv(1-2*P_fa);
    Result = zeros(nBands,1);
    
    for m=1:nBands
        if ((D1(m,1)-D2(m,1)) > Threshold(m,1))
            Result(m)=1;
        else
            Result(m)=0;
        end
    end
    
    nDims=sum(Result);
end

% Purpose: Estimate number of endmembers
% Algorithm name: noise-whitened Harsanyi–Farrand–Chang (HFC) method
% Sensitivity to noise: low
function nDims = NWHFC(imgXYN,t)
    %
    % NWHFC gives the VD number estimated by given false alarm property using
    % NWHFC method.
    %
    % There are two parameters, NWHFC(imgXYN,t) where the imgXYN is the
    % Hyperspectral image cube, which is a 3-D data matrix
    % [XX,YY,bnd] = size(imgXYN), XX YY are the image size,
    % bnd is the band number of the image cube.
    % t is the false alarm probability.
    %
    % HFC uses the NWHFC algorithm developed by Dr. Chein-I Chang,
    % see http://www.umbc.edu/rssipl/. The Matlab code was
    % programmed by Jing Wang in Remote Sensing Signal and
    % Image Processing Lab.
    %

    [XX,YY,bnd] = size(imgXYN);
    pxl_no = XX*YY;
    r = (reshape(imgXYN,pxl_no,bnd))';
    R = (r*r')/pxl_no;
    u = mean(r,2);
    K = R-u*u';
    %======Noise estimation=====
    K_Inverse=inv(K);
    tuta=diag(K_Inverse);
    K_noise=1./tuta;
    K_noise=diag(K_noise);
    %=====Noise whitening===
    y=sqrtm(K_noise) \ r;
    y=reshape(y',XX,YY,bnd);
    %=====Call HFC to estimate===
    nDims=HFC(y,t);

end

% Data Sphering
% Algorithm name: data sphering
% Effectiveness: high
function sphered_data = data_sphering(imgXYN)
    % Initial variables
    [row,column,band]=size(imgXYN);
    % Center the data
    imgXYN_row=reshape(imgXYN,row*column,band);
    imgXYN_zeromean=imgXYN_row-repmat(mean(imgXYN_row),row*column,1);
    cov=imgXYN_zeromean'*imgXYN_zeromean/(row*column);
    % Eigen decomposition
    [V,D]=eig(cov);
    % Transform the data set
    for i=1:band
    sphered_data(i,:)=(V(:,i)'*imgXYN_zeromean')./(D(i,i)^.5*row);
    end
    % Transform the data back
    sphered_data=reshape(sphered_data',row,column,band);

end

% Dimensionality Reduction by Transform
% Algorithm name: principal components analysis (PCA)
% Effectiveness: high
function PCs = PCA(imgXYN,M)
    nBands=size(imgXYN,3);
    nX=size(imgXYN,1);
    nY=size(imgXYN,2);
    x=reshape(imgXYN,nX*nY,nBands);
    x=x';
    L=size(x,1);
    K=size(x,2);
    u=mean(x,2); %dimension of u is 1*L
    x_hat=x-u*ones(1,K);
    m=mean(x,2);
    C=(x*x')/size(x,2)-m*m';
    %===========
    [V,D]=eig(C);
    d=(diag(D))';
    [~,Index]=sort(d);
    for m=1:L
    D_sort(1,m)=d(1,Index(1,Lm));
    V_sort(:,m)=V(:,Index(1,Lm));
    end
    D=diag(D_sort);
    V=V_sort;
    D=D(1:M,1:M);
    V=V(:,1:M);
    %====for the matrix with full column rank, so the
    A=V';
    x_whitened=A*(x_hat);
    PCs=x_whitened;
end

% Dimensionality Reduction by Transform
% Algorithm name: maximum noise fraction (MNF)
% Effectiveness: high
function [ImageCub_MNF,Matrix_of_Vector] = MNF(imgXYN)
    Last_volumn=-10000000000.0000001;
    [height,width,NumberOfSpectrum]=size(imgXYN);
    ImageCub_MNF=zeros(height,width,NumberOfSpectrum);
    %  begin to compute Matrix_of_Vector  %
    meanSpect=zeros(NumberOfSpectrum,1);
    for II=1:height
    for JJ=1:width
    meanSpect=meanSpect+squeeze(imgXYN(II,JJ,:))/height/width;
    end
    end
    TotalCovariance=zeros(NumberOfSpectrum,NumberOfSpectrum);
    for II=1:height
    for JJ=1:width
    TotalCovariance=TotalCovariance+(squeeze(imgXYN(II,JJ,:))-...
    meanSpect)*(squeeze(imgXYN(II,JJ,:))-meanSpect)'/height/width;
    end
    end
    Matrix_F=zeros(NumberOfSpectrum,NumberOfSpectrum);
    Cov_inv=inv(TotalCovariance);
    for II=1:NumberOfSpectrum
    Matrix_F(II,II)=sqrt(Cov_inv(II,II));
    end
    adjusted_Cov=Matrix_F'*TotalCovariance*Matrix_F;
    [V,D]=eig(adjusted_Cov);
    eig_value=zeros(NumberOfSpectrum,2);
    for II=1:NumberOfSpectrum
    eig_value(II,1)=D(II,II);
    eig_value(II,2)=II;
    end
    %disp(eig_value);
    V_sort_min_to_max=sortrows(eig_value,1);
    Matrix_of_Vector_before=zeros(NumberOfSpectrum,NumberOfSpectrum);
    for II=1:NumberOfSpectrum
    Matrix_of_Vector_before(:,II)=squeeze(V(:,V_sort_min_to_max(NumberOf-...
    Spectrum-II+1,2)));
    end
    Matrix_of_Vector=Matrix_F*Matrix_of_Vector_before;
    %  end of computing Matrix_of_Vector  %
    for II=1:height
    for JJ=1:width
    r=squeeze(imgXYN(II,JJ,:))-meanSpect;
    ImageCub_MNF(II,JJ,:)=Matrix_of_Vector'*r;
    end
    end
end

% Dimensionality Reduction by Transform
% Algorithm name: Independent Component Analysis
% Effectiveness: high
function [ICs] = My_fastica_v5(imgXYN,M)
    bnd=size(imgXYN,3);
    xx=size(imgXYN,1);
    yy=size(imgXYN,2);
    x=reshape(imgXYN,xx*yy,bnd);
    x=x';
    L=size(x,1);
    K=size(x,2);
    %====Data sphering =====
    u=mean(x,2); %dimension of u is 1*L
    x_hat=x-u*ones(1,K);
    m=mean(x,2);
    C=(x*x')/size(x,2)-m*m';
    %===========
    [V,D]=eig(C);
    A=inv(sqrtm(D))*V'; % A is the whitening matrix....
    x_whitened=A*(x_hat);
    %=======
    clear x;
    clear x_hat;
    %====rank the eigenvalues, which is used for using the eigenvector as
    %initialization
    %
    % d=(diag(D))';
    % [oo,Index]=sort(d);
    %
    %
    % for m=1:L
    % D_sort(1,m)=d(1,Index(1,L+1-m));
    % V_sort(:,m)=V(:,Index(1,L+1-m));
    % end
    %
    % D=diag(D_sort);
    % V=V_sort;
    %=====
    %====Sphering finished
    threshold = 0.0001;
    B=zeros(L);
    for round=1:M
    fprintf('IC %d', round);
    %===initial condition ===
    w=rand(L,1)-0.5;
    % w=V(:,round); Eigenvectors initialization
    %===
    w=w-B*B'*w;
    w=w/norm(w);
    wOld=zeros(size(w));
    wOld2=zeros(size(w));
    i=1;
    while i<=1000
    w=w-B*B'*w;
    w=w/norm(w);
    fprintf('.');
    if norm(w-wOld)<threshold || norm(w+wOld)<threshold
    fprintf('Convergence after %d steps\n', i);
    B(:,round)=w;
    W(round,:)=w';
    break;
    end
    wOld2=wOld;
    wOld=w;
    w=(x_whitened*((x_whitened'*w).^3))/K-3*w;
    w=w/norm(w);
    i=i+1;
    end
    if (i>1000)
    fprintf('Warning! can not converge after 1000 steps \n, no more components');
    break;
    end
    round=round+1;
    end
    ICs=W*x_whitened;
    % figure;
    % for m=1:M
    % s=reshape(abs(ICs(m,:)),xx,yy); subplot(6,8,m);
    % imagesc(s); axis off;colormap(gray);
    % end

end

% Dimensionality Reduction by Transform
% Algorithm name: Independent Component Analysis SPICA-DR (ICA-DR1)
% Effectiveness: high
function [IC_sorted]=sort_IC_DR1(imgXYN,M)
    clear J;
    bnd=size(imgXYN,3);
    xx=size(imgXYN,1);
    yy=size(imgXYN,2);
    [ICs]=My_fastica_v5(imgXYN,M);
    ICs=abs(ICs);
    % %====Show the ICs in the original order===
    % figure;
    % for m=1:size(ICs,1)
    % s=reshape(abs(ICs(m,:)),xx,yy);
    % s=255*(s-min(min(s))*ones(size(s,1),size(s,2)))/(max(max(s))
%     -min(min(s)));
    % temp=mean(reshape(s,xx*yy,1));
    % subplot(6,8,m); imshow(uint8(s));
    % % title(m);
    % %
    % end
    %======Calculate the contrast function
    for m=1:size(ICs,1)
    s=ICs(m,:);
    var1=var(s);
    mean1=mean(s);
    sita=sqrt(var1);
    skew_temp=sum((s-mean1).^3)/(xx*yy-1);
    kurt_temp=sum((s-mean1).^4)/(xx*yy-1);
    J(1,m)=(skew_temp.^2)/12+((kurt_temp-3).^2)/48;
    end
    %======IC sorting ========
    [~,b]=sort(J);
    b=flipud(b');
    IC_sorted=ICs(b',:);
    %=========Show the ICS after sorting...
    figure;
    for m=1:size(ICs,1)
    s=reshape(abs(IC_sorted(m,:)),xx,yy);
    s=255*(s-min(min(s))*ones(size(s,1),size(s,2)))/(max(max(s)) ...
    -min(min(s)));
    temp=mean(reshape(s,xx*yy,1));
    subplot(6,8,m); imshow(uint8(s));
    % title(m);
    %
    %
    end
end

% Dimensionality Reduction by Transform
% Algorithm name: Independent Component Analysis RICA-DR (ICA-DR2)
% Effectiveness: high
function [IC_selected]=sort_IC_DR2(imgXYN,M,run_times)
    bnd=size(imgXYN,3);
    xx=size(imgXYN,1);
    yy=size(imgXYN,2);
    fprintf('first ICA run \n');
    [ICs]=My_fastica_v5(imgXYN,M*2);
    set1=abs(ICs);
    set_com=set1;
    size1=ones(1);
    for round=1:run_times
    fprintf('ICA run, order is %d \n',round+1);
    [ICs]=My_fastica_v5(imgXYN,2*M);
    set2=abs(ICs);
    set_com_new=[];
    distance=0;
    for m=1:size(set_com,1)
    for n=1:size(set2,1)
    temp1=sqrt(sum(set_com(m,:).^2)); % SAM
    temp2=sqrt(sum(set2(n,:).^2));
    distance(m,n)=acos(sum(set_com(m,:).*set2(n,:))/(temp1*temp2));
    end
    t=distance(m,:)<=0.5;
    if (sum(t)>=1)
    set_com_new=cat(1,set_com_new,set_com(m,:));
    end
    end
    set_com=set_com_new;
    fprintf('the size of set_Com is');
    size(set_com_new)
    size1(round,1)=size(set_com_new,1);
    if (size(set_com_new,1)<=M)
       break;
    end
    end
    IC_selected=set_com;

end

% Dimensionality Reduction by Transform
% Algorithm name: Independent Component Analysis IDICA-DR (ICA-DR3)
% Effectiveness: high
function [Loc,Sig]=My_ATGP(imgXYN,M)
    bnd=size(imgXYN,3);
    xx=size(imgXYN,1);
    yy=size(imgXYN,2);
    r=reshape(imgXYN,xx*yy,bnd);
    r=r';
    %=====Find the first point
    temp=sum(r.*r);
    [~,b]=max(temp);
    if (rem(b,xx)==0)
    Loc(1,1)=b/xx;
    Loc(1,2)=xx;
    elseif (floor(b/xx)==0)
    Loc(1,1)=1;
    Loc(1,2)=b;
    else
    Loc(1,1)=floor(b/xx)+1; % y
    Loc(1,2)=b-xx*floor(b/xx); % x
    end
    Sig(:,1)=r(:,b);
    fprintf('1\n');
    %==========
    for m=2:M
    U=Sig;
    P_U_perl=eye(bnd)-U*((U'*U)\U');
    y=P_U_perl*r;
    temp=sum(y.*y);
    [~,b]=max(temp);
    if (rem(b,xx)==0)
    Loc(m,1)=b/xx;
    Loc(m,2)=xx;
    elseif (floor(b/xx)==0)
    Loc(m,1)=1;
    Loc(m,2)=b;
    else
    Loc(m,1)=floor(b/xx)+1; % y
    Loc(m,2)=b-xx*floor(b/xx); % x
    end
    Sig(:,m)=r(:,b);
    disp(m)
    end
    %
    % figure; imagesc(imgXYN(:,:,30)); colormap(gray); hold on
    % axis off
    % axis equal
    % for m=1:size(Loc,1)
    % plot(Loc(m,1),Loc(m,2),'o','color','g');
    % text(Loc(m,1)+2,Loc(m,2),num2str(m),'color','y','FontSize',12);
    % end
%
end

% Dimensionality Reduction by Transform
% Algorithm name: Independent Component Analysis IDICA-DR (ICA-DR3)
% Effectiveness: high
function ICs = My_fastica_DR3(imgXYN,nComps)
    bnd=size(imgXYN,3);
    xx=size(imgXYN,1);
    yy=size(imgXYN,2);
    x=reshape(imgXYN,xx*yy,bnd);
    x=x';
    L=size(x,1);
    K=size(x,2);
    %====Sphering =====
    u=mean(x,2); %dimension of u is 1*L
    x_hat=x-u*ones(1,K);
    m=mean(x,2);
    C=(x*x')/size(x,2)-m*m';
    %===========
    [V,D]=eig(C);
    A=inv(sqrtm(D))*V'; % A is the whitening matrix....
    x_whitened=A*(x_hat);
    %====for cuprite data===
    clear x;
    clear x_hat;
    %=========for initialization
    [~,Sig]=My_ATGP(reshape(x_whitened',xx,yy,L),nComps);
    W_initial=Sig;
    threshold = 0.0001;
    B=zeros(L);
    %===============find the first point=========
    for round=1:nComps
    fprintf('IC %d', round);
    %===Initial condition
    w=W_initial(:,round);
    %===
    w=w-B*B'*w;
    w=w/norm(w);
    wOld=zeros(size(w));
    wOld2=zeros(size(w));
    i=1;
    while i<=1000
    w=w-B*B'*w;
    w=w/norm(w);
    fprintf('.');
    if norm(w-wOld)<threshold || norm(w+wOld)<threshold
    fprintf('Convergence after %d steps\n', i);
    B(:,round)=w;
    W(round,:)=w';
    break;
    end
    wOld2=wOld;
    wOld=w;
    w=(x_whitened*((x_whitened'*w).^3))/K-3*w;
    w=w/norm(w);
    i=i+1;
    end
    if (i>1000)
    fprintf('Warning! cannot converge after 1000 steps \n, no more components');
    break;
    end
    round=round+1;
    end
    ICs=W*x_whitened;
    figure
    for m=1:nComps
    s=reshape(abs(ICs(m,:)),xx,yy);
    s=255*(s-min(min(s))*ones(size(s,1),size(s,2)))/(max(max(s)) -min(min(s)));
    temp=mean(reshape(s,xx*yy,1));
    subplot(5,6,m); imshow(uint8(s));
    %
    end
    A=W;
    % x_re_ICA=V*sqrtm(D)*(A*ICs)+u*ones(1,xx*yy);

end

% Dimensionality Reduction by Transform
% Algorithm name: high-order statistics DR (HOS-DR)
% Effectiveness: high
function ICs = high_order(imgXYN,nComps,k,Initial)
    % imgXYN is the image cube.
    % M is the number of components to generated using high order
    % k is the order of statistics. For example, k=3, is the skewness, k=4, is
    % the kurtosis, k=5 is the 5th moment, and so on...
    % Initial condition preference, 0 is the random initial, 1 is the eigen
    % initial, 2 is the unity initial. default is 0
    if nargin < 3
    fprintf('Please identify the order of statistics!');
    ICs=[];
    else
    if nargin < 4
    Initial=0;
    end
    bnd=size(imgXYN,3);
    xx=size(imgXYN,1);
    yy=size(imgXYN,2);
    x=reshape(imgXYN,xx*yy,bnd);
    x=x';
    L=size(x,1);
    K=size(x,2);
    %===input x is a matrix with size=L*K;
    %====Sphering =====
    u=mean(x,2); %dimension of u is 1*L
    x_hat=x-u*ones(1,K);
    %===jing's code of cov
    m=mean(x,2);
    C=(x*x')/size(x,2)-m*m';
    %===========
    [V,D]=eig(C);
    A=(sqrtm(D))\V'; % A is the whitening matrix....
    x_whitened=A*(x_hat);
    clear x;
    clear x_hat;
    % Seperating , using high-order
    threshold = 0.01;
    B=zeros(L);
    y=x_whitened;
    W=ones(1,L);
    P_U_perl=eye(bnd);
    for round=1:nComps
        fprintf('IC %d', round);
    %===initial condition ===
    switch (Initial)
    case 0
    w=rand(L,1);
    case 1
    w=V(:,round); % Final version
    case 2
    w=ones(L,1);
    otherwise
    w=rand(L,1);
    end
    i=1;
    while i<=100 % maximum times of trying..
    a=(y.*repmat((w'*y).^(k-2),bnd,1))*y'; %skewness
    a=a/K; % get the sample mean as expectation
    [V,D]=eig(a);
    D=abs(D);
    [~,I]=max(diag(D));
    V1 = V(:,I);
    fprintf('.');
    distance(round,1,i)=norm(w-V1);
    distance(round,2,i)=norm(w+V1);
    if norm(w-V1)<threshold || norm(w+V1)<threshold
    fprintf('Convergence after %d steps\n', i);
    B(:,round)=w;
    W(round,:)=w';
    break;
    end
    w=V1;
    i=i+1;
    end
    %===if not converge. then use the results after 10 iterations
    B(:,round)=w;
    W(round,:)=w';
    %=======
    P_U_perl=eye(bnd)-W'*(W*W')\W;
    y=P_U_perl*y;
    end
    ICs=W*x_whitened;
    figure;
    for m=1:nComps
        s=reshape(abs(ICs(m,:)),xx,yy);
    s=255*(s-min(min(s))*ones(size(s,1),size(s,2)))/(max(max(s)) ...
    -min(min(s)));
    temp=mean(reshape(s,xx*yy,1));
    subplot(5,6,m); imshow(uint8(s));
    %
    end
    end

end

% Endmember Extraction Algorithms
% Algorithm name: Pixel Purity Index
function [eeindex,score,duration] = PPI(imgXYN,nSkewers)
    % The Matlab PPI algorithm
    % ---- Inputs ----
    % 'imagecub' - The hyperspectral image cube
    % 'skewer_no' - The number of skewers
    % ---- Outputs ----
    % 'eeindex' - The locations of the final endmembers (x,y)
    % 'score' - The PPI score of each pixel
    % 'duration - The number of seconds used to run this program
    
    % Initial Variables
    [rows,columns,bands]=size(imgXYN);
    score=zeros(rows*columns,1);
    switch_results=1;
    % Record the start CPU time
    start=cputime();
    % Separate the total number of skewers into several sets and each set uses 500
%     skewers
    skewer_sets=floor(nSkewers/500)+1;
    last_skewer_no=mod(nSkewers,500);
    for i=1:skewer_sets
    if (skewer_sets-i) == 0
    nSkewers=last_skewer_no;
    else
    nSkewers=500;
    end
    % Generate skewers
    rand('state',sum(100*clock));
    skewers=rand(bands,nSkewers)-0.5;
    % Normalize skewers
    for i=1:nSkewers
    skewers(:,i)=skewers(:,i)/norm(skewers(:,i));
    end
    % project every sample vector to the skewers
    projcub=reshape(imgXYN, rows*columns, bands);
    proj_result=projcub*skewers;
    % Find the extrema set for each skewer and add 1 to their score
    for i=1:nSkewers
    max_pos=find(proj_result(:,i)==max(proj_result(:,i)));
    min_pos=find(proj_result(:,i)==min(proj_result(:,i)));
    score(max_pos)=score(max_pos)+1;
    score(min_pos)=score(min_pos)+1;
    end
    end
    % Find the pixel which has score larger than 0
    result=find(score>0);
    % Find the position of the p highest scores
    %result=[];
    %for i=1:skewer_no,
    % result=[result find(max(score)==score,1)];
    % score(find(max(score)==score,1))=0;
    %end
    % Convert one dimension to two dimension index
    if(switch_results)
    eeindex=translate_index(result,rows,columns,1);
    else
    if(mod(result,rows)==0)
    eeindex(2,:)=floor(result./rows);
    else
    eeindex(2,:)=floor(result./rows)+1;
    end
    eeindex(1,:)=mod(result-1,rows)+1;
    end
    duration=cputime()-start;

end

% Endmember Extraction Algorithms
% Algorithm name: Pixel Purity Index
% Effectiveness: high
function [FinalPositions,running_time] = FIPPIoptimized(imgXYN,InitialSkewers)
    % Fast Iterative Pixel Purity Index Algorithm
    %
    % Input parameters:
    % ———————————————
    % Image: Hyperspectral image data after MNF dimensionality reduction
    % InitialSkewers: Positions of ATGP-generated pixels
    %
    % Output parameter:
    % ———————————————
    % FinalPositions: Positions of FIPPI-generated endmember pixels
    % running_time - The total running time used by this run
    %
    % Authors: Chein-I Chang and Antonio Plaza
    % Minor Modified by Chao-Cheng Wu
    % Check CPU time at the beginning
    start=cputime;
    % Code initialization for data and visualization
    [ns,nl]=size(imgXYN);
    Extrema=zeros(ns,nl);
    ProjectionScores=zeros(ns,nl);
    subplot(2,1,1);
    imagesc(imgXYN(:,:,1)); colormap(gray);
    title('Pixels extracted by FPPI:');
    set(gca,'DefaultTextColor','black','xtick',[],'ytick',[],'dataaspectratio',[
    1 1 1]);
    po1 = get(gca,'position');
    % Use ATGP-generated pixels as the initial skewers
    NewSkewers=InitialSkewers;
    % Begin iterative process
    Other = 1;
    SkewersUsed = [];
    while (Other >= 1)
    [ne,~]=size(NewSkewers);
    disp(['Iteration: ' int2str(Other)]);
    disp(['Skewers: ' int2str(ne)]);
    for k = 1:ne
    [ne_old]=size(SkewersUsed);
    skewer=squeeze(imgXYN(NewSkewers(k,1),NewSkewers(k,2),:));
    skewer=skewer/norm(skewer);
    SkewersUsed = union(SkewersUsed,skewer);
    [ne_new]=size(SkewersUsed);
    subplot(2,1,2);
    drawnow;
    plot(skewer);
    title(['Current skewer: ' int2str(k)]);
    if (ne_new~=ne_old)
    % Project all the sample data vectors onto this particular skewer
    for i=1:ns
    for j=1:nl
    pixel = squeeze(imgXYN(i,j,:));
    ProjectionScores(i,j) = dot(skewer,pixel);
    end
    end
    % Obtain the extrema set for each skewer (maximum and minimum
    % projection)
    [~,mpos] = max(ProjectionScores(:));
    [~,pos] = min(ProjectionScores(:));
    mposx = floor((mpos-1)/ns)+1; mposy = mod(mpos-1,ns)+1;
    posx = floor((pos-1)/ns)+1; posy = mod(pos-1,ns)+1;
    % Display the pixel positions of the pixels in the extrema set
    drawnow;
    subplot(2,1,1);
    text(mposx,mposy,'o','Margin',1,'HorizontalAlignment','center',...
    'FontSize',22,'FontWeight','light','FontName','Garamond','Color',...
    'yellow');
    drawnow;
    subplot(2,1,1);
    text(posx,posy,'o','Margin',1,'HorizontalAlignment','center',...
    'FontSize',22,'FontWeight','light','FontName','Garamond','Color',...
    'yellow');
    % Increase PPI count of extrema pixels
    Extrema(posy,posx)=Extrema(posy,posx)+1;
    Extrema(mposy,mposx)=Extrema(mposy,mposx)+1;
    % Incorporate sample vectors with PPI count greater than zero to
    % the skewer set
    vnew = [ mposx mposy ; posx posy ];
    NewSkewers = union(NewSkewers,vnew,'rows');
    end
    end
    % Check stopping rule
    [ne2]=size(NewSkewers);
    if (ne2==ne)
    Other = 0;
    else
    Other = Other+1;
    % Extrema=zeros(ns,nl);
    end
    end
    % Produce the positions of the final endmember set
    Binary = Extrema>0;
    ne = sum(Binary(:));
    disp(['Extracted endmembers: ' int2str(ne)]);
    FinalPositions = zeros(ne,2);
    Current = 1;
    for i=1:ns
    for j=1:nl
    if Binary(i,j)>0
    FinalPositions(Current,1)=i;
    FinalPositions(Current,2)=j;
    Current = Current+1;
    end
    end
    end
    % Check CPU time at the end
    stop=cputime;
    running_time=stop-start;

end

% Endmember Extraction Algorithms
% Algorithm name: N-finder algorithm (N-FINDR)
% Effectiveness: high
function [endmemberindex,duration] = NFINDR(imgXYN,nEndmems)
    % The N-FINDR algorithm
    % ————— Input variables ———————————
    % 'imagecube' - The data transformed components [row column band]
    % 'p' - The number of endmembers to be generated
    %
    % if band > p, then the program will automatically use Singular Value Decomposition
%     to calculate the volume
    % ————— Output variables —————————
    % 'endmemberindex - The locations of the final endmembers (x,y)
    % 'duration - The number of seconds used to run this program
    % Set initial condition
   
    endmemberindex=[];
    newvolume = 0;
    prevolume = -1;
    [row, column, band]=size(imgXYN);
    switch_results=1;
    % Determine to use SVD to calculate the volume or not
    if(band > nEndmems)
    use_svd=1;
    else
    use_svd=0;
    end
    % Start to count the CPU computing time
    start=cputime();
    % Randomly select p initial endmembers
    rand('state',sum(100*clock));
    for i=1:nEndmems
    while(1)
    temp1=round(row*rand);
    temp2=round(column*rand);
    if(temp1>0 && temp2>0)
    break;
    end
    end
    endmemberindex=[endmemberindex;[temp1 temp2]];
    end
    endmemberindex=endmemberindex';
    % Generate endmember vector from reduced cub
    display(endmemberindex);
    endmember=[];
    for i=1:nEndmems
    if(use_svd)
    endmember=[endmember squeeze(imgXYN(endmemberindex(1,i),...
    endmemberindex(2,i),:))];
    else
    endmember=[endmember squeeze(imgXYN(endmemberindex(1,i),...
    endmemberindex(2,i),1:nEndmems-1))];
    end
    end
    % calculate the endmember's volume
    if(use_svd)
    s=svd(endmember);
    endmembervolume=1;
    for i=1:nEndmems
    endmembervolume=endmembervolume*s(i);
    end
    else
    jointmatrix=[ones(1,nEndmems) ; endmember];
    endmembervolume=abs(det(jointmatrix))/factorial(nEndmems-1);
    end
    % The main algorithm
    while newvolume > prevolume % if the new generated endmember volume is larger
%     than the old one, continue the algorithm
    % Use each sample vector to replace the original one, and calculate new volume
    for i=1:row
    for j=1:column
    for k=1:nEndmems
    calc=endmember;
    if(use_svd)
    calc(:,k)=squeeze(imgXYN(i, j, :));
    s=svd(calc);
    volume=1;
    for z=1:nEndmems
    volume=volume*s(z);
    end
    else
    calc(:,k)=squeeze(imgXYN(i, j, 1:nEndmems-1));
    jointmatrix=[ones(1,nEndmems);calculate];
    volume=abs(det(jointmatrix))/factorial(nEndmems-1); % The formula of
    Simplex volume
    end
    if volume > endmembervolume
    endmemberindex(:,k)=[i;j];
    endmember=calc;
    endmembervolume=volume;
    end
    end
    end
    end
    prevolume=newvolume;
    newvolume=endmembervolume;
    end
    stop=cputime();
    duration=stop-start;
    % Switch results for the standard
    if(switch_results)
    endmemberindex(3,:)=endmemberindex(1,:);
    endmemberindex(1,:)=[];
    endmemberindex=endmemberindex';
    end

end

% Endmember Extraction Algorithms
% Algorithm name: simplex growing algorithm (SGA)
% Effectiveness: high
function [endmemberindex,duration] = SGA(imgXYN,nEndmems)
    % Simplex Growing Algorithm
    % - - - - - - - Input variables - - - - - - - - - - - -
    % 'imagecube' - The data transformed components [row column band]
    % 'p' - The number of endmembers to be generated
    %
    % if band > p, then the program will automatically use Singular Value Decomposition
%     to calculate the volume
    % - - - - - - - Output variables - - - - - - - - - - -
    % 'endmemberindex - The locations of the final endmembers (x,y)
    % 'duration - The number of seconds used to run this program
    % Set initial condition
    n=1;
    initial=0;
    [row, column, band]=size(imgXYN);
    % Determine to use SVD to calculate the volume or not
    if(band > nEndmems)
    use_svd=1;
    else
    use_svd=0;
    end
    % Start to count the CPU computing time
    start_time=cputime();
    % Randomly Select a point as the initial point
    endmemberindex=[ceil(row*rand);ceil(column*rand)];
    % The main algorithm
    while n<nEndmems % if get enough endmember group, it stops
    % Generate endmember vector from reduced cub
    endmember=[];
    for i=1:n
    if(use_svd)
    endmember=[endmember squeeze(imgXYN(endmemberindex(1,i),...
    endmemberindex(2,i),:))];
    else
    endmember=[endmember squeeze(imgXYN(endmemberindex(1,i),...
    endmemberindex(2,i),1:n))];
    end
    end
    % Use each sample vector to calculate new volume
    newendmemberindex=[];
    maxvolume=0;
    for i=1:row
    for j=1:column
    if(use_svd)
    jointpoint=[endmember squeeze(imgXYN(i,j,:))];
    s=svd(jointpoint);
    volume=1;
    for z=1:n+1
    volume=volume*s(z);
    end
    else
    jointpoint=[endmember squeeze(imgXYN(i,j,1:n))];
    jointmatrix=[ones(1,n+1);jointpoint];
    volume=abs(det(jointmatrix))/factorial(n); % The formula of a simplex
%     volume
    end
    if volume > maxvolume
    maxvolume=volume;
    newendmemberindex=[i;j];
    end
    end
    end
    endmemberindex=[endmemberindex newendmemberindex]; % Add this pixel into
%     the endmember group
    %nfinder_plot(endmemberindex);
    n=n+1;
    if initial==0 % Use new pixel as the initial pixel
    n=1;
    endmemberindex(:,1)=[];
    initial=initial+1;
    end
    end
    duration=cputime()-start_time;
    % Switch the results back to X and Y
    endmemberindex(3,:)=endmemberindex(1,:);
    endmemberindex(1,:)=[];
    endmemberindex=endmemberindex';

end

% Supervised LSMA and KLSMA
% Algorithm name: orthogonal subspace projection (OSP)
function temp = LSOSP(image,d,U)
    % Least Square Orthogonal Subspace Projection
    % input: image = image cube
    % d = desire signature vector
    % U = undesired signature matrix
    % output: temp = resulting image cube
    [x,y,z]=size(image);
    temp = zeros(x,y);
    %Find the projectors that is orthogonal complement of U
    [l] = size(U);
    I = eye(l,l);
    Pu=I-U*(U'*U)\U';
    lsosp = (d'*Pu)/(d'*Pu*d);
    % perform least-squares-based estimator on all image vectors
    for i = 1:x
    for j = 1:y
    for k = 1:z
    r(k) = image(i,j,k);
    end
    temp(i,j)=lsosp*r';
    end
    end
end

% Supervised LSMA and KLSMA
% Algorithm name: orthogonal subspace projection (OSP)
function temp = KOSP(image,d,U,sig)
    % Kernel based LSOSP function
    % Input:
    % image = image cube input
    % d = desired signature, example: [2;3;4]
    % U = undesired signature matrix
    % sig = parameter that control RBF kernel function
    % output:
    % temp = resulting map
    [x,y,z]=size(image);
    temp = zeros(x,y);
    % perform least squares-based estimator on all image vectors
    KdU = kernelized(d,U,sig,0);%disp(KdU),
    KUU = kernelized(U,U,sig,0);%disp(KUU),
    Kdd = kernelized(d,d,sig,0);
    KUd = kernelized(U,d,sig,0);%disp(KUd),
    for i = 1:x
    for j = 1:y
    for k = 1:z
    r(k,1) = image(i,j,k);
    end
    Kdr = kernelized(d,r,sig,0);%disp(Kdr),
    KUr = kernelized(U,r,sig,0);%disp(KUr),
    temp(i,j)=(Kdr-KdU*(KUU)\KUr);%/(Kdd-KdU*inv(KUU)*KUd);
    end
    end
end

function results = kernelized(x,y,d,chk)
    % kernelization function
    x_l = size(x,2);
    y_l = size(y,2);
    results = zeros(x_l,y_l);
    for i = 1:x_l
    for j = 1:y_l
    results(i,j)= exp((-1/2)*(norm(x(:,i)-y(:,j))^2)/(d^2));
    %RBF kernel (can be changed)
    end
    end
    if chk == 1
    results = results-(sum(sum(results))/(x_l*y_l))*ones(x_l,y_l);
    elseif chk == 2
    N = (1/(x_l*y_l))*ones(x_l,y_l);
    results = results-N*results-results*N+N*results*N;
    end

end

% Supervised LSMA and KLSMA
% Algorithm name: nonnegativity least squares (NCLS)
function [abundance,error_vector] = NCLS(MatrixZ,r1)
    % input MatrixZ is the signatures of endmembers. It is of size [ bands p].
    % input x is the signature whose abundance is to be estimated.
    % output abundance is the abundance of each material in r1. It is of size [p 1].
    % output error_vector is the error vector of size [bands 1].
    % This function is written according to Dr. Chang's first book , P 47
    x=r1; %rename r1 as x;
    M=size(MatrixZ,2);
    count_R=0;
    count_P=M;
    R=zeros(M,1);
    P=ones(M,1);
    %tolerance=0.000001;
    d=zeros(M,1);
    Alpha_ls=(MatrixZ'*MatrixZ)\MatrixZ'*x;
    Alpha_ncls=Alpha_ls;
    min_Alpha_ncls=min(Alpha_ncls);
    M_t_r=MatrixZ'*x;
    invMtM=inv(MatrixZ'*MatrixZ);
    while(min_Alpha_ncls<-0.000000001)
    for II=1:M
    if((Alpha_ncls(II)<0) && (P(II)==1))
    R(II)=1;
    P(II)=0;
    end %%% end of if (Alpha_ncls(II)<0)
    end % end of for II=1:M
    S=R;
    goto_step6=1;
    while(1)
    sum_R=sum(R);
    Alpha_R=zeros(sum_R,1);
    count_for_Alpha_R=0;
    for II=1:M
        if (R(II)==1)
    count_for_Alpha_R=count_for_Alpha_R+1;
    Alpha_R(count_for_Alpha_R)=Alpha_ls(II);
    index_for_Lamda(count_for_Alpha_R)=II;
    end
    end
    count_1_for_P=0;
    Sai_column=[];
    for II=1:M
    if (P(II)~=1)
    Sai_column=[Sai_column squeeze(invMtM(:,II)) ];
    end
    end
    Sai=[];
    for II=1:M
    if (P(II)~=1)
    Sai=[Sai squeeze(Sai_column(II,:)) ];
    end
    end
    Lamda=(Sai)\Alpha_R;
    if(max(Lamda)<0)
    break;
    end
    [~,index_Max_Lamda] = max(Lamda);
    P(index_for_Lamda(index_Max_Lamda))=1;
    R(index_for_Lamda(index_Max_Lamda))=0;
    sum_R=sum(R);
    Alpha_R=zeros(sum_R,1);
    count_for_Alpha_R=0;
    for II=1:M
    if (R(II)==1)
    count_for_Alpha_R=count_for_Alpha_R+1;
    Alpha_R(count_for_Alpha_R)=Alpha_ls(II);
    index_for_Lamda(count_for_Alpha_R)=II;
    end
    end
    Sai_column=[];
    for II=1:M
    if (P(II)~=1)
    Sai_column=[Sai_column squeeze(invMtM(:,II)) ];
    end
    end
    Sai=[];
    for II=1:M
    if (P(II)~=1)
    Sai=[Sai
    squeeze(Sai_column(II,:)) ];
    end
    end
    Lamda=inv(Sai)*Alpha_R;
    Phai_column=[];
    for II=1:M
    if (P(II)~=1)
    Phai_column=[Phai_column squeeze(invMtM(:,II)) ];
    end
    end
    if (size(Phai_column,2)~=0)
    Alpha_s=Alpha_ls-Phai_column*Lamda;
    else
    Alpha_s=Alpha_ls;
    end
    goto_step6=0;
    find_smallest_in_S=zeros(M,2);
    find_smallest_in_S(:,1)=Alpha_s;
    find_smallest_in_S(:,2)=[1:M]';
    sort_find=sortrows(find_smallest_in_S,1);
    for II=1:M
    if ((S(II)==1)&&(Alpha_s(II)<0))
    P(II)=0;
    R(II)=1;
    goto_step6=1;
    end
    end
    end % end of while (gotostep6==1)
    Phai_column=[];
    for II=1:M
    if (P(II)~=1)
    Phai_column=[Phai_column squeeze(invMtM(:,II)) ];
    end
    end
    if (size(Phai_column,2)~=0)
    Alpha_ncls=Alpha_ls-Phai_column*Lamda;
    else
    Alpha_ncls=Alpha_ls;
    end
    min_Alpha_ncls=min(Alpha_ncls);
    end % end of while
    abundance=zeros(M,1);
    for II=1:M
    if (Alpha_ncls(II)>0)
    abundance(II)=Alpha_ncls(II);
    end
    end
    error_vector=MatrixZ*abundance-x;

end

