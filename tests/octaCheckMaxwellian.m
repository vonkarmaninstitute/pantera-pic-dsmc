close all
clear
clc

% Load Maxwellian
dd = load('../OUTPMAXW');

vx = dd(:,1);
vy = dd(:,2);
vz = dd(:,3);

[N,X] = hist(vx,200);

figure
plot(X, N, 'k', 'linewidth', 2)
hold on
plot(-X, N, 'r', 'linewidth', 2)

