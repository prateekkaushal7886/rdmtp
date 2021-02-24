vall=zeros(20,1);
valx=zeros(20,1);

for k = 1:1:20
   %% CREATE/LOAD GEOMETRY
% Flag to specify creating or loading airfoil
flagAirfoil.XFoilCreate = 1;                                                
flagAirfoil.XFoilLoad   = 0;                                                


Vinf = 1;                                                                   
AoA  = 0;  
nac=2*1000+4*100+(k-1)*5;
NACA =int2str(nac) ;                                                              
alpha = AoA*(pi/180);                                                       

% Plotting flags
flagPlot = [0;          % Airfoil with panel normal vectors
            0;          % Geometry boundary pts, control pts, first panel, second panel
            0;          % Cp vectors at airfoil surface panels
            1;          % Pressure coefficient comparison (XFOIL vs. VPM)
            0;          % Airfoil streamlines
            0];         % Pressure coefficient contour

%% XFOIL - CREATE/LOAD AIRFOIL

% PPAR menu options
PPAR.N  = '160';                                                            % "Number of panel nodes"
PPAR.P  = '4';                                                              % "Panel bunching parameter"
PPAR.T  = '1.5';                                                            % "TE/LE panel density ratios"
PPAR.R  = '1';                                                              % "Refined area/LE panel density ratio"
PPAR.XT = '1 1';                                                            % "Top side refined area x/c limits"
PPAR.XB = '1 1';                                                            % "Bottom side refined area x/c limits"

% Call XFOIL function to obtain the following:
% - Airfoil coordinates
% - Pressure coefficient along airfoil surface
% - Lift, drag, and moment coefficients
[xFoilResults,success] = XFOIL(NACA,PPAR,AoA,flagAirfoil);                  % Get XFOIL results for prescribed airfoil
if (success == 0)                                                           % If user canceled airfoil dialog box
    return;                                                                 % Exit the program
end

% Separate out results from XFOIL function results
afName  = xFoilResults.afName;                                              % Airfoil name
xFoilX  = xFoilResults.X;                                                   % X-coordinate for Cp result
xFoilY  = xFoilResults.Y;                                                   % Y-coordinate for Cp result
xFoilCP = xFoilResults.CP;                                                  % Pressure coefficient
XB      = xFoilResults.XB;                                                  % Boundary point X-coordinate
YB      = xFoilResults.YB;                                                  % Boundary point Y-coordinate
xFoilCL = xFoilResults.CL;                                                  % Lift coefficient
xFoilCD = xFoilResults.CD;                                                  % Drag coefficient
xFoilCM = xFoilResults.CM;                                                  % Moment coefficient


numPts = length(XB);                                                        
numPan = numPts - 1;         

% Check for direction of points

edge = zeros(numPan,1);                                                     
for i = 1:1:numPan                                                          
    edge(i) = (XB(i+1)-XB(i))*(YB(i+1)+YB(i));                              
end
sumEdge = sum(edge);                                                        
fprintf(NACA);
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
    if (phiD(i) < 0)                                                        % Make all panel angles positive [deg]
        phiD(i) = phiD(i) + 360;
    end
      
end

% Compute angle of panel normal w.r.t horizontal and include AoA
deltaD             = phiD + 90;                                             % Angle of panel normal
betaD              = deltaD - AoA;                                          % Angle of panel normal and AoA [deg]
betaD(betaD > 360) = betaD(betaD > 360) - 360;                              % Make sure angles aren't greater than 360 [deg]


phi  = phiD.*(pi/180);                                                      
beta = betaD.*(pi/180);                                                     

% Geometric integral (normal [K] and tangential [L])

[K,L] = COMPUTE_KL_VPM(XC,YC,XB,YB,phi,S); 





% Populate A matrix
A = zeros(numPan,numPan);                                                   
for i = 1:1:numPan                                                          
    for j = 1:1:numPan                                                      
        if (j == i)                                                         
            A(i,j) = 0;                                                     
        else                                                                
            A(i,j) = -K(i,j);                                               
        end
    end
end

% Populate b array

b = zeros(numPan,1);                                                        
for i = 1:1:numPan                                                          
    b(i) = -Vinf*2*pi*cos(beta(i));                                         
end

% Satisfy the Kutta condition

A(numPan,:)   = 0;                                                          
A(numPan,1)   = 1;                                                          
A(numPan,end) = 1;                                                          
b(numPan)     = 0;

% Compute gamma values from system of equations
gamma = A\b;                                    

%% COMPUTE PANEL VELOCITIES AND PRESSURE COEFFICIENTS

% Compute velocities on each panel
Vt = zeros(numPan,1);                                                       
Cp = zeros(numPan,1);                                                       
for i = 1:1:numPan                                                          
    addVal  = 0;                                                            
    for j = 1:1:numPan                                                      
        addVal = addVal - (gamma(j)/(2*pi))*L(i,j);                         
    end
    
    Vt(i) = Vinf*sin(beta(i)) + addVal + gamma(i)/2;                        
    Cp(i) = 1-(Vt(i)/Vinf)^2;                                               
end

%% COMPUTE LIFT AND MOMENT COEFFICIENTS

% Compute normal and axial force coefficients
CN = -Cp.*S.*sin(beta);                                                     
CA = -Cp.*S.*cos(beta);                                                     

% Compute lift and drag coefficients
CL = sum(CN.*cosd(AoA)) - sum(CA.*sind(AoA));                               
CM = sum(Cp.*(XC-0.25).*S.*cos(phi));                                       

%% PLOTTING

T = linspace(0,360,1000);                                                  % Angle array to compute dashed circle
x = cosd(T);                                                                % Circle X points
y = sind(T);                                                                % Circle Y points

% Plot shape polygon with panel normal vectors
figure(k);                                                                 % Create figure
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
zoom reset;         


vall(k)=CL;
valx(k)=xFoilCL;

end;
    