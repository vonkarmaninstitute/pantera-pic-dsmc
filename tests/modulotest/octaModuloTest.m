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
hold on
plot([xMin, xMax, xMax, xMin, xMin], [yMin, yMin, yMax, yMax, yMin], 'k', 'linewidth', 2)
axis equal
grid on
hold on
plot(xP, yP, 'or', 'linewidth',3)
plot(xP1, yP1, 'om', 'linewidth',3)


