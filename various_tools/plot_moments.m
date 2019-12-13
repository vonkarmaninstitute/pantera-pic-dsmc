close all
clear
clc

% Load moments file
dd = load('../dumps/global_moments.dat');

tt  = dd(:,1);
rho = dd(:,2);
ux  = dd(:,3);
uy  = dd(:,4);
uz  = dd(:,5);
Pxx = dd(:,6);
Pxy = dd(:,7);
Pxz = dd(:,8);
Pyy = dd(:,9);
Pyz = dd(:,10);
Pzz = dd(:,11);
qx  = dd(:,12);
qy  = dd(:,13);
qz  = dd(:,14);


figure
plot(tt, ux, 'r', 'linewidth', 2)
hold on
plot(tt, uy, 'g', 'linewidth', 2)
plot(tt, uz, 'b', 'linewidth', 2)
xlabel('Time [s]')
ylabel('Velocity [m/s]')

figure
plot(tt, Pxx, 'c', 'linewidth', 2)
hold on
plot(tt, Pxy, 'm', 'linewidth', 2)
plot(tt, Pxz, 'r', 'linewidth', 2)
plot(tt, Pyy, 'g', 'linewidth', 2)
plot(tt, Pyz, 'b', 'linewidth', 2)
plot(tt, Pzz, 'r', 'linewidth', 2)
xlabel('Time [s]')
ylabel('Pressure tensor [Pa]')

figure
plot(tt, qx, 'r', 'linewidth', 2)
plot(tt, qy, 'g', 'linewidth', 2)
plot(tt, qz, 'b', 'linewidth', 2)
xlabel('Time [s]')
ylabel('Heat flux vector [W/m2 I guess]')

n = rho/9.109e-31;

T = (Pxx + Pyy + Pzz)./n/3/1.38e-23;

figure
plot(tt, T*1.38e-23/1.602e-19, 'b', 'linewidth', 2)
ylabel('T [eV]')
xlabel('Time [s]')
title('Memento excitatio...')
pbaspect([2,1,1])

