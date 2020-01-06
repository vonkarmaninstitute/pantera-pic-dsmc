close all
clear
clc

% Fields
E = 10000; % [V/m]
B = 0.01; % [T]

Evec = [E; 0; 0];
Bvec = [0; 0; B];

q = -1.602e-19; % [C]
M = 9.1e-31;   % [kg]

omega_c = q/M*B;

% 
Nt = 800;
tt = linspace(0, 8*(2*pi/abs(omega_c)), Nt)';
dt = tt(2) - tt(1); 

x0 = 0;
y0 = 0;

vx0 = 1e6;
vy0 = -2e6;

% Init vectors - numerical solution
xx_euler = zeros(numel(tt),1);
yy_euler = zeros(numel(tt),1);

vx_euler = zeros(numel(tt),1);
vy_euler = zeros(numel(tt),1);

xx_euler(1) = x0;
yy_euler(1) = y0;

vx_euler(1) = vx0;
vy_euler(1) = vy0;

% Init vectors - numerical solution - Boris
xx_Boris = zeros(numel(tt),1);
yy_Boris = zeros(numel(tt),1);

vx_Boris = zeros(numel(tt),1);
vy_Boris = zeros(numel(tt),1);
vz_Boris = zeros(numel(tt),1);

xx_Boris(1) = x0;
yy_Boris(1) = y0;

vx_Boris(1) = vx0 + q/M*(E + vy0*B)*dt/2; % Advect Boris velocity by dt/2
vy_Boris(1) = vy0 - q/M*vx0*B*dt/2;       % Advect Boris velocity by dt/2

% Init vectors - analytical solution
xx_analyt = zeros(numel(tt),1);
yy_analyt = zeros(numel(tt),1);

vx_analyt = zeros(numel(tt),1);
vy_analyt = zeros(numel(tt),1);

xx_analyt(1) = x0;
yy_analyt(1) = y0;

% Init velocity with half timestep
vx_analyt(1) = vx0;
vy_analyt(1) = vy0;

% Start time loop
for ii = 1:numel(tt)-1

  % ++++++++++ Forward Euler ++++++++++
  vx_euler(ii+1) = vx_euler(ii) + dt*q/M*(E + vy_euler(ii)*B);
  vy_euler(ii+1) = vy_euler(ii) - dt*q/M*(vx_euler(ii)*B);

  xx_euler(ii+1) = xx_euler(ii) + vx_euler(ii+1)*dt;
  yy_euler(ii+1) = yy_euler(ii) + vy_euler(ii+1)*dt;

  % ++++++++++ Analytical ++++++++++
  psi = omega_c*dt;
  vx_analyt(ii+1) = -sin(psi)*(-vy_analyt(ii) - E/B) + cos(psi)*vx_analyt(ii);
  vy_analyt(ii+1) = -cos(psi)*(-vy_analyt(ii) - E/B) - sin(psi)*vx_analyt(ii) - E/B;

  xx_analyt(ii+1) = xx_analyt(ii) + 1/omega_c*(-(1-cos(psi))*(-vy_analyt(ii) - E/B) + sin(psi)*vx_analyt(ii));
  yy_analyt(ii+1) = yy_analyt(ii) - 1/omega_c*(sin(psi)*(-vy_analyt(ii) - E/B) + (1-cos(psi))*vx_analyt(ii)) - psi/omega_c*E/B;

  % +++++++++++ Boris method ++++++++++
  % Advect position
  xx_Boris(ii+1) = xx_Boris(ii) + vx_Boris(ii)*dt;
  yy_Boris(ii+1) = yy_Boris(ii) + vy_Boris(ii)*dt;

  % Advect velocity
  v_ii = [vx_Boris(ii); vy_Boris(ii); vz_Boris(ii)];
  v_minus = v_ii + q*Evec/M*dt/2;
  tvec    = q*Bvec/M*dt/2;
  v_prime = v_minus + cross(v_minus, tvec);
  svec    = 2*tvec/(1+dot(tvec,tvec));
  v_plus  = v_minus + cross(v_prime, svec);

  vv_iip1 = v_plus + q*Evec/M*dt/2;

  vx_Boris(ii+1) = vv_iip1(1);
  vy_Boris(ii+1) = vv_iip1(2);
  vz_Boris(ii+1) = vv_iip1(3);

end

% Plot
figure
subplot(1,2,1)
plot(vx_euler, vy_euler, 'b', 'linewidth', 2)
hold on
plot(vx_analyt, vy_analyt, '-or', 'linewidth', 2)
plot(vx_Boris, vy_Boris, '--k', 'linewidth', 2)
xlabel('vx')
ylabel('vy')

subplot(1,2,2)
plot(xx_euler, yy_euler, 'b', 'linewidth', 2)
hold on
plot(xx_analyt, yy_analyt, '-or', 'linewidth', 2)
plot(xx_Boris, yy_Boris, '--k', 'linewidth', 2)
xlabel('x')
ylabel('y')
legend('Forward Euler', 'Analytical', 'Boris', "location", 'south')
title(['Dt = ', num2str(dt*abs(omega_c)), ' 1/ \omega_c'])
