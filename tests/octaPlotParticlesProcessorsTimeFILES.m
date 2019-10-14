close all
clear
clc

page_screen_output(0);

% Parameters
timesteps = [1:1:100];
Nproc = 9;

% Domain box
XX = [0, 2.1];
YY = [-0.5,0.5];
ZZ = [-0.1,0.1];

colors = 'rgbkmyc';

figure

for tid = timesteps(1):timesteps(end)

  hold off
  for iproc = 0:Nproc-1
 
    filename = sprintf('../dumps/proc_%05d_time_%08d', iproc, tid);

    % File may be empty
    try
      dd = load(filename);
    catch 
      dd = [NaN, NaN, NaN, NaN, NaN];
    end

    col_now = colors(mod(iproc,7)+1);
    plot3(dd(:,2), dd(:,3), dd(:,4), '.', 'color', col_now);
    % plot3(dd(1,2), dd(1,3), dd(1,4), '.', 'color', col_now); % Only plot first particle
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
  zlabel('z [m]')

%   imgname = sprintf('./IMGS/img_%05d.png',tid); 
%   print(imgname)
  pause(0.1)

end

