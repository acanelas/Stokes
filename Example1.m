% Analysis of Example 1

EST.PasoMax = 0.01;
EST.TOL = 1e-10; % Tolerance
EST.MaxIt = 400; % Maximum number of iterations
EST.MinIt = 3;   % Minimum number of iterations

% functions
EST.fun = @funFU;   % Unmodified version
EST.gfun = @gfunFU; % Unmodified version derivative

% Remeshing parameters
EST.REM = false; % remeshing disabled in this example
EST.Delta = 0.05;
EST.Beta  = 2*sin(pi/40); % 40 nodes in a circle

% Geometry
NE1 = 60; % initial number of elements of interior boundary
NE2 = 80; % initial number of elements of exterior boundary
Rs = 0.4; % initial radius of interior boundary
Rg = 0.8; % initial radius of exterior boundary

% Gravity parameter b
EST.bA = -1.0; % Vertical magnitude

% Velocity g
EST.gA = 4.031531359547362e-02; % Vertical magnitude
EST.gB = 0.0;                   % Rotational magnitude

% Pressure parameter
EST.p0 = 1.0;

% Material parameters
EST.Gamma = 0.0; % Parameter gamma
EST.Mu = 1.0;    % Parameter mu

% Initial mesh
t1 = (0:-2*pi/NE1:-(2*pi-1e-5))';

% Coord X;  Coord Y;
EST.MatNod1 = [
    Rs*cos(t1),  Rs*sin(t1);
    ];

% Nodes;
EST.MatEle1 = [
    (1:NE1)', (2:NE1+1)';
    ];

t2 = (0:2*pi/NE2:(2*pi-1e-5))';
% Coord X;  Coord Y;
EST.MatNod2 = [
    Rg*cos(t2),  Rg*sin(t2);
    ];

% Nodes;
EST.MatEle2 = [
    (1:NE2)', (2:NE2+1)';
    ];

EST.MatEle1(NE1,2) = 1;
EST.MatEle2(NE2,2) = 1;

% Initial u
EST.XU = zeros(2*NE2,1);

% Initial p
EST.XP = zeros(2*NE1,1);

% Function g
EST.XG = zeros(2*NE1,1);
EST.XG(1:2:end) = 0      - EST.gB*Rs*sin(t1);
EST.XG(2:2:end) = EST.gA + EST.gB*Rs*cos(t1);

% Initial guess x0 --------------------------------------------------------
NumVar = 2*NE1 + 3*NE2 + 1;
EST.x0 = zeros(NumVar,1);
EST.x0(1:2*NE1)             = EST.XP;
EST.x0(2*NE1+1:2*NE1+2*NE2) = EST.XU;

% Analysis ----------------------------------------------------------------
figurita(EST);
EST = NewAnalysis(EST);
figurita(EST);
