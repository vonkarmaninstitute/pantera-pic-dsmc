close all
clear
clc

% Domain
xMin = -0.9221;
xMax = 0.823;

yMin = 0.1;
yMax = 2.2;

Lx = xMax - xMin;
Ly = yMax - yMin;

% Particle
xP = -1.1;
yP = 0.5;

xP1 = xMin + mod(xP - xMax, Lx);
yP1 = yP;



figure

for tID = 1:1000
  plot([xMin, xMax, xMax, xMin, xMin], [yMin, yMin, yMax, yMax, yMin], 'k', 'linewidth', 2)
  axis equal
  hold on
  grid on
  
  plot(xP1, yP1, 'om', 'linewidth',3)
  pause(0.05)
  hold off

  xP1 = xP1 - 0.03;
  yP1 = yP1 - 0.02;

  xP1 = xMin + mod(xP1 - xMax, Lx);
  yP1 = yMin + mod(yP1 - yMax, Ly);
end




