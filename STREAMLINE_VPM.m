function [Nx,Ny] = STREAMLINE_VPM(XP,YP,XB,YB,phi,S)


% Number of panels
numPan = length(XB)-1;                                                      % Number of panels (control points)

% Initialize arrays
Nx = zeros(numPan,1);                                                       % Initialize Nx integral array
Ny = zeros(numPan,1);                                                       % Initialize Ny integral array

% Compute Nx and Ny
for j = 1:1:numPan                                                          % Loop over all panels
    % Compute intermediate values
    A  = -(XP-XB(j))*cos(phi(j))-(YP-YB(j))*sin(phi(j));                    % A term
    B  = (XP-XB(j))^2+(YP-YB(j))^2;                                         % B term
    Cx = sin(phi(j));                                                       % Cx term (X-direction)
    Dx = -(YP-YB(j));                                                       % Dx term (X-direction)
    Cy = -cos(phi(j));                                                      % Cy term (Y-direction)
    Dy = XP-XB(j);                                                          % Dy term (Y-direction)
    E  = sqrt(B-A^2);                                                       % E term
    if (~isreal(E))
        E = 0;
    end
    
    % Compute Nx
    term1 = 0.5*Cx*log((S(j)^2+2*A*S(j)+B)/B);                              % First term in Nx equation
    term2 = ((Dx-A*Cx)/E)*(atan2((S(j)+A),E) - atan2(A,E));                 % Second term in Nx equation
    Nx(j) = term1 + term2;                                                  % Compute Nx integral
    
    % Compute Ny
    term1 = 0.5*Cy*log((S(j)^2+2*A*S(j)+B)/B);                              % First term in Ny equation
    term2 = ((Dy-A*Cy)/E)*(atan2((S(j)+A),E) - atan2(A,E));                 % Second term in Ny equation
    Ny(j) = term1 + term2;                                                  % Compute Ny integral
    
	% Zero out any NANs, INFs, or imaginary numbers
    if (isnan(Nx(j)) || isinf(Nx(j)) || ~isreal(Nx(j)))
        Nx(j) = 0;
    end
    if (isnan(Ny(j)) || isinf(Ny(j)) || ~isreal(Ny(j)))
        Ny(j) = 0;
    end
end
