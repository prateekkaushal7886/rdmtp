
%Input Parameters
naca=0015;
TE=1;                       % 1 for open trailing edge and 2 for closed trailing edge
gridPoints=100;
gridtype=1;                     %1 for uniform and 2 for non-uniform



m=floor(naca/1000);
p=rem(floor((naca/100)),10);
t=rem(naca,100);


M = m/100;
P = p/10;
T = t/100;

a0 = 0.2969;
a1 = -0.1260;
a2 = -0.3516;
a3 = 0.2843;
if (TE == 1)
    a4 = -0.1015;                                                        % Open trailing edge
elseif (TE == 2)
    a4 = -0.1036;                                                           % Closed trailing edge
end


%X coordinate points
if (gridtype == 1)                                                       % Uniform spacing
    x = linspace(0,1,gridPoints)';
elseif (gridtype == 2)                                                   % Non-uniform spacing
    beta = linspace(0,pi,gridPoints)';
    x = (0.5*(1-cos(beta)));
end

% Camber line and camber line gradient
yc     = ones(gridPoints,1);
dyc_dx = ones(gridPoints,1);
theta  = ones(gridPoints,1);
for i = 1:1:gridPoints
    if (x(i) >= 0 && x(i) < P)
        yc(i)     = (M/P^2)*((2*P*x(i))-x(i)^2);
        dyc_dx(i) = ((2*M)/(P^2))*(P-x(i));
    elseif (x(i) >=P && x(i) <=1)
        yc(i)     = (M/(1-P)^2)*(1-(2*P)+(2*P*x(i))-(x(i)^2));
        dyc_dx(i) = ((2*M)/((1-P)^2))*(P-x(i));
    end
    theta(i) = atan(dyc_dx(i));
end


% Thickness distribution
yt = 5*T.*((a0.*sqrt(x)) + (a1.*x) + (a2.*x.^2) + (a3.*x.^3) + (a4.*x.^4));

% Upper surface points
xu = x(:)  - yt(:).*sin(theta);
yu = yc(:) + yt(:).*cos(theta);

% Lower surface points
xl = x(:) + yt(:).*sin(theta);
yl = yc(:) - yt(:).*cos(theta);

% Camber line X points
xc = x;


figure(1);
        cla; hold on; grid on;
        set(gcf,'Color','White');
         plot(xu,yu,'k-'); 
          plot(xl,yl,'k-');
          xlim([0 1]);
          %legend([p1,p2],{'Data 1','Data 2'});                                      % Add legend
xlabel('X Units');                                                          % Set X-label
ylabel('Y Units');

% Set Y-label
axis equal;                                                                 % Set axes equal
zoom reset;         


