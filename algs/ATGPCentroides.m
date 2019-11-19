function [U,P] = ATGPCentroides(centroides,p)
% ATGP algorithm for endmember extraction.
% ------------------------------------------------------------------------------
% Input:   HIM : hyperspectral image data in matrix form (pixels by bands).
%                Pixels by rows.
%
% Output:  U   : set of extracted endmembers [nchannels x p] 
%          P   : Number of the extracted endmembers.
% 
% Copyright (2007) GRNPS group @ University of Extremadura, Spain. 

%%disp(' === Start ATGP run ===')

% Obtener tiempo CPU actual | get current CPU time
t1=cputime;

% Obtener tama�o de la imagen (muestras,lineas,bandas) | get image size
[nc,nb]=size(centroides);


% Visualizar imagen | visualize image
%%imagesc(mean(HIM,3)); colormap(gray); 
%%set(gca,'DefaultTextColor','black','xtick',[],'ytick',[],'dataaspectratio',[1 1 1])
%%po1 = get(gca,'position');

% Calculo del pixel (vector) con mayor intensidad en la imagen | 
% Calculate the pixel (vector) with major intensity in the image


max = 0;
pos =0;
for i = 1:nc
       r = squeeze( centroides(i,:) );
       bright = r'*r;
       if bright > max
           max = bright;
           pos = i;
       end
end

centro = centroides(pos,:);

% max = 0;
% centro(nb) = 0;
% for i = 1:nc
%     centro(:) = centro(:) + reshape(centroides(i,:),nb,1 ) ;
% end
% centro(:) = centro(:) / nc;

%%max = -999999;
%%%pos = 0;
%%for i = 1:nc
%%    angulo = dot(squeeze(centroides(i,:)), centro(:)) / (norm(squeeze(centroides(i,:)),2 )*norm(centro(:),2) ) ;
%%    if (angulo > max);
%%        max = angulo;
 %%       pos = i;
 %%   end
%%end




% El pixel con mas intensidad es el pixel inicial del proceso |
% The pixel with more intensity is the initial pixel of the process
t0 = squeeze(centro(:));

% Calculo de la matriz identidad | Generate the identity matrix.
I = eye(nb,nb);

% Inicializacion de la matriz de pixels puros |
% Initialization of the pure pixels matrix
U = [];
U = [U t0];

% Inicializacion de la matriz de posiciones |
% Initialization of the positions matrix
P = zeros(p);
%%P(1)=pos;

%%disp(sprintf('. found pixel @ coordinates x=%5d & y=%5d',posx,posy))

% Visualizacion de la posicion del primer pixel seleccionado |
% Visualization of the position of the first chosen pixel
%%drawnow;
%text(posy,posx,'o','Margin',1,'HorizontalAlignment','center','FontSize',22,'FontWeight','light','FontName','Garamond','Color','yellow');
%%text(posy,posx,'o','Color','yellow');

% Algoritmo ATGP
P(1,1) = pos;
for i = 2:p
    UC = U(:,1:i-1);
    % Calculo de la proyeccion ortogonal con respecto a los pixels
    % actualmente seleccionados. Esta parte puede sustituirse por cualquier
    % otra distancia |
    % Calculate the orthogonal projection with respect to the pixels at present chosen. 
    % This part can be replaced with any other distance
    PU = I-UC*pinv(UC'*UC)*UC';
    maximum = 0;
    % Calculo del pixel mas distinto a los ya seleccionados en funcion de la
    % proyeccion ortogonal (o cualquier otra distancia que se seleccione) |
    % Calculate the pixel most different from the already selected ones according to
    % the orthogonal projection (or any other distance selected)
    for n = 1:nc
        if (centroides(n,1) ~= 0 )
            r = squeeze( centroides(n,:) );
            result = PU*r';
            val = result'*result;
            if (val > maximum)
                maximum = val;
                pos = n
            end
        end
    end
    % El siguiente pixel seleccionado es el mas diferente a los ya seleccionados |
    % The next chosen pixel is the most different from the already chosen ones
    ti = squeeze( centroides(pos,:) )';
    % Mostrar posiciones de dicho pixel por pantalla |
    % Show positions of the above mentioned pixel at screen
    %%disp(sprintf('. found pixel @ coordinates x=%5d & y=%5d',posx,posy))
    % Almacenar posiciones en matriz de posiciones |
    % Store positions in the matrix of positions
    P(i,1)=pos;
    % Almacenar pixel en matriz de pixels puros |
    % Store positions in the matrix of pure pixels
    U = [U ti];
    % Visualizar pixel seleccionado en pantalla |
    % Visualize pixel selected on screen
    %%drawnow;
    %% text(pos,'o','Margin',1,'HorizontalAlignment','center','FontSize',22,'FontWeight','light','FontName','Garamond','Color','yellow');
    %%text(posy,posx,'o','Color','yellow');

end

% Obtener tiempo CPU actual | get current CPU time
U = U(:,2:p);
t2=cputime;

% Mostrar tiempo total en ejecucion del algoritmo |
% Show total execution time of the algorithm
%% disp( sprintf('. Total CPU processing time .................... %6.3f [s]  ', (t2-t1) ) );
%%disp(' === Eng ATGP ===');
end
