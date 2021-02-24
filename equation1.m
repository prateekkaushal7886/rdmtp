
clear;
clc;

%% CREATE/LOAD GEOMETRY
                                                          
load = 'Airfoil';                                                            
AoA  = 0;
Vinf = 1; 

fileName = 'naca2415';                                                    
[XB,YB]  = LOAD_AIRFOIL_SELIG(fileName);



% Number of panels
numPan = length(XB)-2; 
XB(numPan+2)=[];
YB(numPan+2)=[];

%% CHECK PANEL DIRECTIONS - FLIP IF NECESSARY

% Check for direction of points
edge = zeros(numPan,1);                                                     % Initialize edge value array
for i = 1:1:numPan                                                          % Loop over all panels
    edge(i) = (XB(i+1)-XB(i))*(YB(i+1)+YB(i));                              % Compute edge values
end
sumEdge = sum(edge);                                                        % Sum all edge values

% If panels are CCW, flip them (don't if CW)
if (sumEdge < 0)                                                            % If panels are CCW
    XB = flipud(XB);                                                        % Flip the X-data array
    YB = flipud(YB);                                                        % Flip the Y-data array
end

%% PANEL METHOD GEOMETRY - REF [1]

% Initialize variables
XC   = zeros(numPan,1);                                                     % Initialize control point X-coordinate array
YC   = zeros(numPan,1);                                                     % Initialize control point Y-coordinate array
S    = zeros(numPan,1);                                                     % Initialize panel length array
phiD = zeros(numPan,1);                                                     % Initialize panel orientation angle array [deg]

% Find geometric quantities of the airfoil
for i = 1:1:numPan                                                          % Loop over all panels
    XC(i)   = 0.5*(XB(i)+XB(i+1));                                          % X-value of control point
    YC(i)   = 0.5*(YB(i)+YB(i+1));                                          % Y-value of control point
    dx      = XB(i+1)-XB(i);                                                % Change in X between boundary points
    dy      = YB(i+1)-YB(i);                                                % Change in Y between boundary points
    S(i)    = (dx^2 + dy^2)^0.5;                                            % Length of the panel
	phiD(i) = atan2d(dy,dx);                                                % Angle of the panel (positive X-axis to inside face) [deg]
    if (phiD(i) < 0)                                                        % Make all panel angles positive [deg]
        phiD(i) = phiD(i) + 360;
    end
end

% Compute angle of panel normal w.r.t horizontal and include AoA
deltaD             = phiD + 90;                                             % Angle from positive X-axis to outward normal vector [deg]
betaD              = deltaD - AoA;                                          % Angle between freestream vector and outward normal vector [deg]
betaD(betaD > 360) = betaD(betaD > 360) - 360;                              % Make all panel angles between 0 and 360 [deg]

% Convert angles from [deg] to [rad]
phi  = phiD.*(pi/180);                                                      % Convert from [deg] to [rad]
beta = betaD.*(pi/180);                                                     % Convert from [deg] to [rad]

%% COMPUTE VORTEX PANEL STRENGTHS - REF [5]

% Geometric integral (normal [K] and tangential [L])
% - Refs [2] and [3]
[K,L] = COMPUTE_KL_VPM(XC,YC,XB,YB,phi,S);                                  % Compute geometric integrals

% Populate A matrix
% - Simpler option: A = -K;
A = zeros(numPan,numPan);                                                   % Initialize the A matrix
for i = 1:1:numPan                                                          % Loop over all i panels
    for j = 1:1:numPan                                                      % Loop over all j panels
        if (j == i)                                                         % If the panels are the same
            A(i,j) = 0;                                                     % Set A equal to zero
        else                                                                % If panels are not the same
            A(i,j) = -K(i,j);                                               % Set A equal to negative geometric integral
        end
    end
end

% Populate b array
% - Simpler option: b = -Vinf*2*pi*cos(beta);
b = zeros(numPan,1);                                                        % Initialize the b array
for i = 1:1:numPan                                                          % Loop over all panels
    b(i) = -Vinf*2*pi*cos(beta(i));                                         % Compute RHS array
end

% Satisfy the Kutta condition
pct    = 100;                                                               % Panel replacement percentage
panRep = floor((pct/100)*numPan);                                           % Replace this panel with Kutta condition eqn.
if (panRep == 0)                                                            % If we specify the first panel
    panRep = 1;                                                             % Make sure the index is not zero
end
A(panRep,:)   = 0;                                                          % Set all columns of the replaced panel equal to zero
A(panRep,1)   = 1;                                                          % Set first column of replaced panel equal to 1
A(panRep,end) = 1;                                                          % Set last column of replaced panel equal to 1
b(panRep)     = 0;                                                          % Set replaced panel value in b array equal to zero

% Compute gamma values from system of equations
gamma = A\b;                                                                % Compute all vortex strength values

%% COMPUTE PANEL VELOCITIES AND PRESSURE COEFFICIENTS

% Compute velocities on each panel
Vt = zeros(numPan,1);                                                       % Initialize tangential velocity array
Cp = zeros(numPan,1);                                                       % Initialize pressure coefficient array
for i = 1:1:numPan                                                          % Loop over all i panels
    addVal  = 0;                                                            % Reset the summation value to zero
    for j = 1:1:numPan                                                      % Loop over all j panels
        addVal = addVal - (gamma(j)/(2*pi))*L(i,j);                         % Sum all tangential vortex panel terms
    end
    
    Vt(i) = Vinf*sin(beta(i)) + addVal + gamma(i)/2;                        % Compute tangential velocity by adding uniform flow and i=j terms
    Cp(i) = 1-(Vt(i)/Vinf)^2;                                               % Compute pressure coefficient
end

%% COMPUTE LIFT AND MOMENT COEFFICIENTS

% Compute normal and axial force coefficients
CN = -Cp.*S.*sin(beta);                                                     % Normal force coefficient []
CA = -Cp.*S.*cos(beta);                                                     % Axial force coefficient []

% Compute lift and drag coefficients
CL = sum(CN.*cosd(AoA)) - sum(CA.*sind(AoA));                               % Decompose axial and normal to lift coefficient []
CM = sum(Cp.*(XC-0.25).*S.*cos(phi));                                       % Moment coefficient []

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
zoom reset;         




