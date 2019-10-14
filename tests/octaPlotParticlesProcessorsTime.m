close all
clear
clc

page_screen_output(0);

% Domain box
XX = [0, 2.1];
YY = [-0.5,0.5];
ZZ = [-0.1,0.1];

% Load particles
%dd = load('outp2');
%dd = load('../outp');
%dd = load('outt');
dd = load('../outt');

Nproc = max(dd(:,5)) + 1;
Ntime = max(dd(:,1)) + 1;

colors = 'rgbkmycrgbkmycrgbkmycrgbkmycrgbkmycrgbkmyc'

figure

for tid = 1:Ntime

  hold off  

  ID_time_now = find(dd(:,1) == tid);
  
  dd_now = dd(ID_time_now, :);

  for iproc = 0:Nproc-1
  
    col_now = colors(iproc+1);
  
    IDplotNow = find(dd_now(:,5) == iproc);
  
    plot3(dd_now(IDplotNow, 2), dd_now(IDplotNow, 3), dd_now(IDplotNow, 4), '.', 'color', col_now)
    hold on
  
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

  pause(0.1)

end


