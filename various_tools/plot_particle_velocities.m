close all
clear
clc

page_screen_output(0);

files_list = dir('../dumps/proc_00000_time_*');
% files_list = dir('/media/starlight/Maxtor/PANTERA_data/proc_00000_time_*');

dt = 5e-8;    % Timestep

for ii = 1:numel(files_list)

  % Load file
  dd = load(['../dumps/',files_list(ii).name]);
  fprintf('Data from: %s\n', files_list(ii).name)

  t_vec(ii)  = dd(4,1);
  vx_vec(ii) = dd(4,5);
  vy_vec(ii) = dd(4,6);
  vz_vec(ii) = dd(4,7);

end 

figure 
plot(t_vec, vy_vec, '-+b', 'linewidth', 2)

figure
plot(vy_vec, vz_vec, 'or', 'linewidth', 2)

figure
hold off
for ii = 1:10:numel(files_list)

  hold off
  plot(vy_vec(1:ii), vz_vec(1:ii), 'r', 'linewidth', 2)
  hold on
  plot(vy_vec(ii), vz_vec(ii), 'or', 'linewidth', 2)
  pause(0.1)
end
