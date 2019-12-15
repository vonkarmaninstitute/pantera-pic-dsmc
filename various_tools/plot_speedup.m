close all
clear
clc

tt = [21.4, 13.5, 10.1, 7.4, 3.6, 2.9, 2.3]; % times [s]
Np = [1,    2,    3,    4,   8,   12,  16];  % Number of proc

figure
plot(Np, 1./(tt/tt(1)), '-ob', 'linewidth', 2)
xlabel('Number of processors')
ylabel('Speedup')
