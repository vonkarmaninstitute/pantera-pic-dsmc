function Ndot = MaxwellianFlux(n, T, RGAS, u)

beta = 1/sqrt(2*RGAS*T);

s = u*beta;

Ndot = n/beta*(exp(-s^2) + sqrt(pi)*s*(1 + erf(s)))/(2*sqrt(pi))

return
