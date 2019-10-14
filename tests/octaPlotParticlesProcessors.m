close all
clear
clc

% Domain box
XX = [0, 2.1];
YY = [-0.5,0.5];
ZZ = [-0.1,0.1];

% Load particles
%dd = load('outp2');
%dd = load('../outp');
%dd = load('outt');
dd = load('../outpp');

Nproc = max(dd(:,4)) + 1;

colors = 'rgbkmycrgbkmycrgbkmycrgbkmycrgbkmycrgbkmyc'

figure
hold on
for iproc = 0:Nproc-1

  col_now = colors(iproc+1)

  IDplotNow = find(dd(:,4) == iproc);

  plot3(dd(IDplotNow, 1), dd(IDplotNow, 2), dd(IDplotNow, 3), '.', 'color', col_now)

end

axis equal

plot3([XX(1), XX(2), XX(2), XX(1), XX(1)], [YY(1), YY(1), YY(2), YY(2), YY(1)], [ZZ(1), ZZ(1), ZZ(1), ZZ(1), ZZ(1)], 'k')
plot3([XX(1), XX(2), XX(2), XX(1), XX(1)], [YY(1), YY(1), YY(2), YY(2), YY(1)], [ZZ(2), ZZ(2), ZZ(2), ZZ(2), ZZ(2)], 'k')
plot3([XX(1), XX(1)], [YY(1), YY(1)], [ZZ(1), ZZ(2)], 'k')
plot3([XX(2), XX(2)], [YY(1), YY(1)], [ZZ(1), ZZ(2)], 'k')
plot3([XX(1), XX(1)], [YY(2), YY(2)], [ZZ(1), ZZ(2)], 'k')
plot3([XX(2), XX(2)], [YY(2), YY(2)], [ZZ(1), ZZ(2)], 'k')

xlabel('x [m]')
ylabel('y [m]')
