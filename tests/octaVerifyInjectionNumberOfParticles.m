close all
clear
clc

% Load simulation result
dd = load('../outp');

% Geometrical parameters
ninject = 1.0e23; % [1/m^3]
uinject = 100.0;  % [m/s]
Tinject = 500.0;  % [m/s]
RGAS    = 200.0;  % [J kg/K ??]

%Ainject = 1*0.2;  % [m^2]  ! Injection along Y (side X = 0)
Ainject = 2.1*0.2;  % [m^2] ! Injection along X (side Y = 0)

vol     = 2.1*1*0.2; % [m^3]
Fnum    = 1.0e18;    % [-]
dt      = 1.e-4;     % [s]
Nt      = numel(dd); % Number of timesteps

tvect = [0:dt:(Nt-1)*dt];

% Compute number of real particles
Nreal = dd*Fnum;

% Compute theoretical flux of particles through surface from Maxwellian VDF
Ndot = Ainject*MaxwellianFlux(ninject, Tinject, RGAS, uinject);

figure
plot(tvect, Nreal,'-or', 'linewidth', 2)
hold on
plot(tvect, Ndot*tvect,'-b', 'linewidth', 2)
xlabel('Timestep')
ylabel('Number of (real) particles in domain')
legend('Nsimulated', 'Ntheoretical')
 
