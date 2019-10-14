close all
clear
clc

dd = load('particles.dat');

figure
plot3(dd(:,1), dd(:,2), dd(:,3), 'or', 'linewidth', 2)
