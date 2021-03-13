
clear;
clc;

%% KNOWNS

% Flag to specify creating or loading airfoil
flagAirfoil.XFoilCreate = 1;                                                
flagAirfoil.XFoilLoad   = 0;                                                


Vinf = 1;                                                                   
AoA  = 0;                                                                   
NACA = '2418';                                                              
alpha = AoA*(pi/180);                                                       

% Plotting flags
flagPlot = [1;          % Airfoil with panel normal vectors
            1;          % Geometry boundary pts, control pts, first panel, second panel
            1;          % Cp vectors at airfoil surface panels
            1;          % Pressure coefficient comparison (XFOIL vs. VPM)
            1;          % Airfoil streamlines
            1];         % Pressure coefficient contour

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

%% CHECK PANEL DIRECTIONS - FLIP IF NECESSARY

% Check for direction of points
edge = zeros(numPan,1);                                                     
for i = 1:1:numPan                                                          
    edge(i) = (XB(i+1)-XB(i))*(YB(i+1)+YB(i));                              
end
sumEdge = sum(edge);                                                        
% If panels are CCW, flip them (don't if CW)
if (sumEdge < 0)                                                            
    XB = flipud(XB);                                                        
    YB = flipud(YB);                                                        
end

%% PANEL METHOD GEOMETRY - REF [1]

% Initialize variables
XC   = zeros(numPan,1);                                                    
YC   = zeros(numPan,1);                                                    
Cp_x=zeros(numPan,1);
S    = zeros(numPan,1);                                                    
phiD = zeros(numPan,1);                                                    

% Find geometric quantities of the airfoil
for i = 1:1:numPan                                                         
    XC(i)   = 0.5*(XB(i)+XB(i+1));                                         
    YC(i)   = 0.5*(YB(i)+YB(i+1));                                         
    dx      = XB(i+1)-XB(i);                                               
    dy      = YB(i+1)-YB(i);                                               
    S(i)    = (dx^2 + dy^2)^0.5;                                           
	phiD(i) = atan2d(dy,dx);                                               
    if (phiD(i) < 0)                                                       
        phiD(i) = phiD(i) + 360;
    end
end


deltaD             = phiD + 90;                                             
betaD              = deltaD - AoA;                                          
betaD(betaD > 360) = betaD(betaD > 360) - 360;                              

phi  = phiD.*(pi/180);                                                      
beta = betaD.*(pi/180);                                                     

%% COMPUTE VORTEX PANEL STRENGTHS - REF [5]


[K,L] = COMPUTE_KL_VPM(XC,YC,XB,YB,phi,S);                                  % Compute geometric integrals


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


b = zeros(numPan,1);                                                        
for i = 1:1:numPan                                                          
    b(i) = -Vinf*2*pi*cos(beta(i));                                         
end

% Satisfy the Kutta condition
pct    = 100;                                                               
panRep = floor((pct/100)*numPan);                                           
if (panRep == 0)                                                            
    panRep = 1;                                                             
end
A(panRep,:)   = 0;                                                          
A(panRep,1)   = 1;                                                          
A(panRep,end) = 1;                                                          
b(panRep)     = 0;                                                          

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
CN = -Cp.*S.*sin(beta);                                                     % Normal force coefficient []
CA = -Cp.*S.*cos(beta);                                                     % Axial force coefficient []

% Compute lift and drag coefficients
CL = sum(CN.*cosd(AoA)) - sum(CA.*sind(AoA));                               % Decompose axial and normal to lift coefficient []
CM = sum(Cp.*(XC-0.25).*S.*cos(phi));                                       % Moment coefficient []


%% COMPUTE STREAMLINES - REF [4]

if (flagPlot(5) == 1 || flagPlot(6) == 1)                                   % If we are plotting 5 or 6
    % Grid parameters
    nGridX = 150;                                                           % X-grid for streamlines and contours
    nGridY = 150;                                                           % Y-grid for streamlines and contours
	xVals  = [-0.5; 1.5];                                                   % X-grid extents [min, max]
    yVals  = [-0.3; 0.3];                                                   % Y-grid extents [min, max]
    
    % Streamline parameters
    stepsize = 0.01;                                                        % Step size for streamline propagation
    maxVert  = nGridX*nGridY*100;                                           % Maximum vertices
    slPct    = 25;                                                          % Percentage of streamlines of the grid
    Ysl      = linspace(yVals(1),yVals(2),floor((slPct/100)*nGridY))';      % Create array of Y streamline starting points
    
    % Generate the grid points
    Xgrid   = linspace(xVals(1),xVals(2),nGridX)';                          % X-values in evenly spaced grid
    Ygrid   = linspace(yVals(1),yVals(2),nGridY)';                          % Y-values in evenly spaced grid
    [XX,YY] = meshgrid(Xgrid,Ygrid);                                        % Create meshgrid from X and Y grid arrays
    
    % Initialize velocities
    Vx = zeros(nGridX,nGridY);                                              % Initialize X-velocity matrix
    Vy = zeros(nGridX,nGridY);                                              % Initialize Y-velocity matrix
    
    % Solve for grid point X and Y velocities
    for m = 1:1:nGridX
        for n = 1:1:nGridY
            XP = XX(m,n);                                                   % Current iteration's X grid point
            YP = YY(m,n);                                                   % Current iteration's Y grid point
            [Nx,Ny] = STREAMLINE_VPM(XP,YP,XB,YB,phi,S);                    % Compute Nx and Ny geometric integrals
            
            [in,on] = inpolygon(XP,YP,XB,YB);
            if (in == 1 || on == 1)                                         % If the grid point is in or on the airfoil
                Vx(m,n) = 0;                                                % Set X-velocity equal to zero
                Vy(m,n) = 0;                                                % Set Y-velocity equal to zero
            else                                                            % If the grid point is outside the airfoil
                Vx(m,n) = Vinf*cosd(AoA) + sum(-gamma.*Nx./(2*pi));         % Compute X-velocity
                Vy(m,n) = Vinf*sind(AoA) + sum(-gamma.*Ny./(2*pi));         % Compute Y-velocity
            end
        end
    end
    
    % Compute grid point velocity magnitude and pressure coefficient
    Vxy  = sqrt(Vx.^2 + Vy.^2);                                             % Compute magnitude of velocity vector []
    CpXY = 1-(Vxy./Vinf).^2;                                                % Pressure coefficient []
end



%% PLOTTING

% FIGURE: Airfoil with panel normal vectors
if (flagPlot(1) == 1)
    figure(1);                                                              % Create figure
    cla; hold on; grid off;                                                 % Get ready for plotting
    set(gcf,'Color','White');                                               % Set color to white
    set(gca,'FontSize',12);                                                 % Set font size
    fill(XB,YB,'k');                                                        % Plot airfoil
    for i = 1:1:numPan                                                      % Loop over all panels
        X(1) = XC(i);                                                       % Set X start of panel orientation vector
        X(2) = XC(i) + S(i)*cosd(betaD(i)+AoA);                             % Set X end of panel orientation vector
        Y(1) = YC(i);                                                       % Set Y start of panel orientation vector
        Y(2) = YC(i) + S(i)*sind(betaD(i)+AoA);                             % Set Y end of panel orientation vector
        plot(X,Y,'r-','LineWidth',2);                                       % Plot panel normal vector
    end
    xlabel('X Units');                                                      % Set X-label
    ylabel('Y Units');                                                      % Set Y-label
	xlim('auto');                                                           % Set X-axis limits to auto
    ylim('auto');                                                           % Set Y-axis limits to auto
    axis equal;                                                             % Set axes equal
    zoom reset;                                                             % Reset zoom
end



% FIGURE: Cp vectors at airfoil control points
if (flagPlot(3) == 1)
    figure(3);                                                              % Create figure
    cla; hold on; grid on;                                                  % Get ready for plotting
    set(gcf,'Color','White');                                               % Set color to white
    set(gca,'FontSize',12);                                                 % Set font size
    Cps = abs(Cp*0.25);                                                     % Scale and make positive all Cp values
    for i = 1:1:length(Cps)                                                 % Loop over all panels
        X(1) = XC(i);                                                       % Control point X-coordinate
        X(2) = XC(i) + Cps(i)*cosd(betaD(i)+AoA);                           % Ending X-value based on Cp magnitude
        Y(1) = YC(i);                                                       % Control point Y-coordinate
        Y(2) = YC(i) + Cps(i)*sind(betaD(i)+AoA);                           % Ending Y-value based on Cp magnitude
        
        if (Cp(i) < 0)                                                      % If pressure coefficient is negative
            p{1} = plot(X,Y,'r-','LineWidth',2);                                   % Plot as a red line
        elseif (Cp(i) >= 0)                                                 % If pressure coefficient is zero or positive
            p{2} = plot(X,Y,'b-','LineWidth',2);                                   % Plot as a blue line
        end
    end
    fill(XB,YB,'k');                                                        % Plot the airfoil as black polygon
	legend([p{1},p{2}],{'Negative Cp','Positive Cp'});                      % Show legend
    xlabel('X Units');                                                      % Set X-label
    ylabel('Y Units');                                                      % Set Y-label
    xlim('auto');                                                           % Set X-axis limits to auto
    ylim('auto');                                                           % Set Y-axis limits to auto
    axis equal;                                                             % Set axes equal
    zoom reset;                                                             % Reset zoom
end

% FIGURE: Pressure coefficient comparison (XFOIL vs. VPM)
if (flagPlot(4) == 1)
    figure(4);                                                              % Create figure
    cla; hold on; grid on;                                                  % Get ready for plotting
    set(gcf,'Color','White');                                               % Set color to white
    set(gca,'FontSize',12);                                                 % Set font size
    midIndX = floor(length(xFoilCP)/2);                                     % Airfoil middle index for XFOIL data
    midIndS = floor(length(Cp)/2);                                          % Airfoil middle index for SPM data
    pXu = plot(xFoilX(1:midIndX),xFoilCP(1:midIndX),'b-','LineWidth',2);    % Plot Cp for upper surface of airfoil from XFOIL
    pXl = plot(xFoilX(midIndX+1:end),xFoilCP(midIndX+1:end),'r-',...        % Plot Cp for lower surface of airfoil from XFOIL
                    'LineWidth',2);
    pVl = plot(XC(1:midIndS),Cp(1:midIndS),'b--o','MarkerFaceColor','r');     % Plot Cp for upper surface of airfoil from SPM
    pVu = plot(XC(midIndS+1:end),Cp(midIndS+1:end),'r--o',...                 % Plot Cp for lower surface of airfoil from SPM
                    'MarkerFaceColor','b');
    legend([pXu,pXl,pVu,pVl],...                                            % Show legend
           {'XFOIL Upper','XFOIL Lower','VPM Upper','VPM Lower'});
    xlabel('X Coordinate');                                                 % Set X-label
    ylabel('Cp');                                                           % Set Y-label
    xlim([0 1]);                                                            % Set X-axis limits
    ylim('auto');                                                           % Set Y-axis limits to auto
    set(gca,'Ydir','reverse')                                               % Reverse direction of Y-axis
    title(['Airfoil: ' xFoilResults.afName ...                              % Title
           ', CL_{VPM}/CL_{XFOIL} = ' ...
           num2str(CL,4) '/' num2str(xFoilCL,4)]);
    zoom reset;                                                             % Reset zoom
end

% FIGURE: Airfoil streamlines
if (flagPlot(5) == 1)
    figure(5);                                                              % Create figure
    cla; hold on; grid on;                                                  % Get ready for plotting
    set(gcf,'Color','White');                                               % Set color to white
    set(gca,'FontSize',12);                                                 % Set font size
    for i = 1:1:length(Ysl)                                                 % Loop over all Y streamline starting points
        sl = streamline(XX,YY,Vx,Vy,xVals(1),Ysl(i),[stepsize,maxVert]);    % Plot the streamline
        set(sl,'LineWidth',2);                                              % Set streamline line width
    end
    fill(XB,YB,'k');                                                        % Plot airfoil as black polygon
	xlabel('X Units');                                                      % Set X-label
    ylabel('Y Units');                                                      % Set Y-label
    xlim(xVals);                                                            % Set X-axis limits
	axis equal;                                                             % Set axes equal
    ylim(yVals);                                                            % Set Y-axis limits
    zoom reset;                                                             % Reset zoom
end

% FIGURE: Pressure coefficient contour
if (flagPlot(6) == 1)
    figure(6);                                                              % Create figure
    cla; hold on; grid on;                                                  % Get ready for plotting
    set(gcf,'Color','White');                                               % Set color to white
    set(gca,'FontSize',12);                                                 % Set font size
    contourf(XX,YY,CpXY,100,'EdgeColor','none');                            % Plot Cp contour
    fill(XB,YB,'k');                                                        % Plot airfoil as black polygon
    xlabel('X Units');                                                      % Set X-label
    ylabel('Y Units');                                                      % Set Y-label
    xlim(xVals);                                                            % Set X-axis limits
	axis equal;                                                             % Set axes equal
    ylim(yVals);                                                            % Set Y-axis limits
    zoom reset;                                                             % Reset zoom
end
