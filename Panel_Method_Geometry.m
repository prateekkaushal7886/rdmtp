
clear;
clc;

%% CREATE/LOAD GEOMETRY
                                                          
load = 'Airfoil';                                                            
AoA  = 0;                                                                   

fileName = 'naca2408';                                                    
[XB,YB]  = LOAD_AIRFOIL_SELIG(fileName);


% Number of panels
numPan = length(XB)-1;                                                      

% Check for direction of points
edge = zeros(numPan,1);                                                     
for i = 1:1:numPan                                                          
    edge(i) = (XB(i+1)-XB(i))*(YB(i+1)+YB(i));                              
end
sumEdge = sum(edge);                                                        

% If panels are CCW, flip them (don't if CW)
if (sumEdge < 0)                                                            
    fprintf('Points are counter-clockwise.  Flipping.\n');                  
    XB = flipud(XB);                                                        
    YB = flipud(YB);                                                        
elseif (sumEdge > 0)                                                        
    fprintf('Points are clockwise.  Not flipping.\n');                      
end

%% COMPUTE GEOMETRIC VARIABLES

% Find geometric quantities of airfoil
XC   = zeros(numPan,1);                                                     
YC   = zeros(numPan,1);                                                     
S    = zeros(numPan,1);                                                     
phiD = zeros(numPan,1);                                                     
for i = 1:1:numPan                                                         
    XC(i)   = 0.5*(XB(i)+XB(i+1));                                          
    YC(i)   = 0.5*(YB(i)+YB(i+1));                                          
    dx      = XB(i+1)-XB(i);                                                
    dy      = YB(i+1)-YB(i);                                                
    S(i)    = (dx^2 + dy^2)^0.5;                                            
	phiD(i) = atan2d(dy,dx);                                                
    
end

% Compute angle of panel normal w.r.t horizontal and include AoA
deltaD             = phiD + 90;                                             % Angle of panel normal
betaD              = deltaD - AoA;                                          % Angle of panel normal and AoA [deg]
betaD(betaD > 360) = betaD(betaD > 360) - 360;                              % Make sure angles aren't greater than 360 [deg]

%% PLOTTING

T = linspace(0,360,1000);                                                  % Angle array to compute dashed circle
x = cosd(T);                                                                % Circle X points
y = sind(T);                                                                % Circle Y points

% Plot shape polygon with panel normal vectors
figure(21);                                                                 % Create figure
cla; hold on; grid off;                                                     % Get ready for plotting
set(gcf,'Color','White');                                                   % Set color to white
set(gca,'FontSize',12);                                                     % Set font size
fill(XB,YB,'k');                                                            % Plot polygon

for i = 1:1:numPan                                                          % Loop over all panels
    X(1) = XC(i);                                                           % Panel starting X point
    X(2) = XC(i) + S(i)*cosd(deltaD(i));                                    % Panel ending X point
    Y(1) = YC(i);                                                           % Panel starting Y point
    Y(2) = YC(i) + S(i)*sind(deltaD(i));                                    % Panel ending Y point
    if (i == 1)                                                             % For first panel
        p1 = plot(X,Y,'r-','LineWidth',2);                                  % Plot first panel normal vector
    elseif (i == 2)                                                         % For second panel
        p2 = plot(X,Y,'r-','LineWidth',2);                                  % Plot second panel normal vector
    else                                                                    % For every other panel
        plot(X,Y,'r-','LineWidth',2);                                       % Plot panel normal vector
    end
end
legend([p1,p2],{'Data 1','Data 2'});                                      % Add legend
xlabel('X Units');                                                          % Set X-label
ylabel('Y Units');                                                          % Set Y-label
axis equal;                                                                 % Set axes equal
zoom reset;                                                                 % Reset zoom


