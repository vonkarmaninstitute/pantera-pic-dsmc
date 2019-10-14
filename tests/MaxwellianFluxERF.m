close all
clear
clc

n = 1.0e22;
T = 500.0;
RGAS = 200.0;

u = 100.0;

beta = 1/sqrt(2*RGAS*T);

s = u*beta;

Ndot = n/beta*(exp(-s^2) + sqrt(pi)*s*(1 + erf(s)))/(2*sqrt(pi))


