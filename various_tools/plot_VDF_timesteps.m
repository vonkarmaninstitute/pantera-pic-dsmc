close all
clear
clc

page_screen_output(0);

files_list = dir('../dumps/proc_00000_time_*');
% files_list = dir('/media/starlight/Maxtor/PANTERA_data/proc_00000_time_*');

dt = 5e-8;    % Timestep

PLOTVDF = figure();

for ii = 1:numel(files_list)

  % Load file
  dd = load(['../dumps/',files_list(ii).name]);
  fprintf('Data from: %s\n', files_list(ii).name)

%   t_vec(ii)  = dd(4,1);
%   vx_vec(ii) = dd(4,5);
%   vy_vec(ii) = dd(4,6);
%   vz_vec(ii) = dd(4,7);

  vy = dd(:,6); 
  vz = dd(:,7);

  % Compute histogram
  [ff, vyP, vzP] = hist2d([vy, vz], 50, 50);

  figure(PLOTVDF)
  surf(vyP, vzP, ff)
  xlim([0, 8e7])
  ylim([0, 8e7])
  pause(0.1)
end 
