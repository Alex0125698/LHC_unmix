%%% Uncomment to run example

% [grid,x,y] = simple2Dexample();
% [X,Y] = meshgrid(x,y);
% pcolor(X,Y,grid(:,:)');

% [grid,x,y,z] = simple3Dexample();
% [X,Y] = meshgrid(x,y);
% y2 = yields(grid);
% subplot(2,2,1);
% surf(X,Y,grid(:,:,10)');
% subplot(2,2,2);
% plot(x,y2{1},'-',x,y2{1},'o');
% subplot(2,2,3);
% plot(y,y2{2},'-',y,y2{2},'o');
% subplot(2,2,4);
% plot(z,y2{3},'-',z,y2{3},'o');

[grid,x,y] = random2Dexample();
[X,Y] = meshgrid(x,y);
y2 = yields(grid);
subplot(2,2,1);
surf(X,Y,grid(:,:)');
subplot(2,2,2);
plot(x,y2{1},'-',x,y2{1},'o');
subplot(2,2,3);
plot(x,y2{2},'-',x,y2{2},'o');

function [grid,x,y] = simple2Dexample

    nDims = 2;
    VdimMins = zeros(nDims,1);
    VdimMaxs = ones(nDims,1);
    VdimPoints = [50;30];
    Vscale = [0.5;0.1];
    Mrotate = [0,pi/8;0,0];
    Vtrans = [-0.3;-0.5];

    grid = gaussianDataGen(nDims,VdimMins,VdimMaxs,VdimPoints,Vscale,Mrotate,Vtrans);
    x = (0:VdimPoints(1)-1)/(VdimPoints(1)-1)*VdimMaxs(1) + VdimMins(1);
    y = (0:VdimPoints(2)-1)/(VdimPoints(2)-1)*VdimMaxs(2) + VdimMins(2);

end

function [grid,x,y,z] = simple3Dexample

    nDims = 3;
    VdimMins = zeros(nDims,1);
    VdimMaxs = ones(nDims,1);
    VdimPoints = [50;30;20];

    S(1).strength = 0.9;
    S(1).Vscale = [0.65;0.15;0.1];
    S(1).Mrotate = [0,pi/4,pi/6;pi/40,0,0;0,0,0];
    S(1).Vtrans = [-0.8;-0.6;-0.40];

    S(2).strength = 1.1;
    S(2).Vscale = [0.3;0.6;0.1];
    S(2).Mrotate = [0,0,0;0,0,0;0,0,0];
    S(2).Vtrans = [-0.20;-0.20;-0.5];

    grid = MultiGaussianDataGen(nDims,VdimMins,VdimMaxs,VdimPoints,S);
    x = (0:VdimPoints(1)-1)/(VdimPoints(1)-1)*VdimMaxs(1) + VdimMins(1);
    y = (0:VdimPoints(2)-1)/(VdimPoints(2)-1)*VdimMaxs(2) + VdimMins(2);
    z = (0:VdimPoints(3)-1)/(VdimPoints(3)-1)*VdimMaxs(3) + VdimMins(3);
    
end

function [grid,x,y] = random2Dexample
    
    nGaussians = 6; % number of gaussians to generate
    nDims = 2;
    VdimMins = zeros(nDims,1);
    VdimMaxs = ones(nDims,1);
    VdimPoints = 60*ones(nDims,1);
  
    for i=1:nGaussians
        S(i).strength = rand(1)+0.5;
        S(i).Vscale = rand(nDims,1)*0.3 + 0.05;
        S(i).Mrotate = rand(nDims)*2*pi;
        S(i).Vtrans = -(0.8*rand(nDims,1)+0.1);   
    end

	grid = MultiGaussianDataGen(nDims,VdimMins,VdimMaxs,VdimPoints,S);    
    x = (0:VdimPoints(1)-1)/(VdimPoints(1)-1)*VdimMaxs(1) + VdimMins(1);
    y = (0:VdimPoints(2)-1)/(VdimPoints(2)-1)*VdimMaxs(2) + VdimMins(2);
        
end

%%% Yields

function y = yields(grid)

    dimLengths = size(grid);
    nDims = length(dimLengths);
    
    for i=1:nDims
        % sum over all dimensions other than i
        y{i} = reshape(sum(grid,[1:i-1;i+1:nDims]),dimLengths(i),1);        
    end

end

%%% Gaussian Generators

function grid = MultiGaussianDataGen(nDims,VdimMins,VdimMaxs,VdimPoints,S)

    grid = zeros(VdimPoints');
    
    for i=1:length(S)
        tmp = gaussianDataGen(nDims,VdimMins,VdimMaxs ...
               ,VdimPoints, S(i).Vscale, S(i).Mrotate, S(i).Vtrans) ;
        grid = grid + S(i).strength*tmp;
    end
    
end

function grid = gaussianDataGen(nDims,VdimMins,VdimMaxs,VdimPoints,Vscale,Mrotate,Vtrans)
    % generate a multidimensional grid of data where each gridpoint
    % has accociated number of yeilds. Gaussian yeild
    % distributions are added to the grid according to the input
    % parameters
    
    % number of roataion matrices
    nRotMats = nDims*(nDims-1)/2;
    % the roataion matrices
    Mrots = eye(nDims+1) + zeros(nDims+1,nDims+1,nRotMats);
    
    counter = 1;
    for row = 1:nDims
        for col = row+1:nDims
            Mrots(row,row,counter) = cos(Mrotate(row,col));
            Mrots(col,col,counter) = cos(Mrotate(row,col));
            Mrots(row,col,counter) = sin(Mrotate(row,col));
            Mrots(col,row,counter) = -sin(Mrotate(row,col));
            counter = counter + 1;
        end
    end

    % Translation matrix
    Mtrans = eye(nDims+1);
    Mtrans(1:nDims,nDims+1) = Vtrans;

    % scaling matrix
    Mscale = diag([1./Vscale;1]);

    % generate vector of all points in desired grid
    nPoints = prod(VdimPoints);
    gridPoints = ones(nDims+1,nPoints);
    
    for row = 1:nDims % coord of each dim
        elemReps = prod([1;VdimPoints((1:row-1)')]);
        rowReps = nPoints / (elemReps*VdimPoints(row));        
        gridPoints(row,:) = repmat(repelem((0:VdimPoints(row)-1),1,elemReps),1,rowReps);
    end
    
    gridPoints = gridPoints .* ([VdimMaxs;1] ./ [(VdimPoints-1);1]) + [VdimMins;0];
    
    % reverse transformations to get into original coordinates
    Mnet = Mscale;   
    for i=1:nRotMats
        Mnet = Mnet * Mrots(:,:,i);
    end 
    Mnet = Mnet * Mtrans;
    
    origGridPoints = Mnet * gridPoints;
    
    % gaussian generator; input vector of arbitrary dimentionality
    gauss = @(x) exp(-sum(x.*x));
    
    % calculate the gaussian at each original coordinate
    gaussWeights = gauss(origGridPoints);
    
    grid = reshape(gaussWeights,VdimPoints');    

end
