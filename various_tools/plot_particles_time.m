close all
clear
clc

page_screen_output(0);

% files_list = dir('../dumps/proc_00000_time_*');
% files_list = dir('/media/starlight/Maxtor/PANTERA_data/proc_00000_time_*');
files_list = dir('/media/starlight/Maxtor/PANTERA_data_test_Boris/proc_00000_time_*');

dt = 1e-11;    % Timestep

figure
for ii = 1:numel(files_list)

  % Load file
  % dd = load(['../dumps/',files_list(ii).name]);
  dd = load(['/media/starlight/Maxtor/PANTERA_data_test_Boris/',files_list(ii).name]);
  fprintf('Data from: %s\n', files_list(ii).name)

  t_vec(ii)  = dd(1,1);
  vx_vec(ii) = dd(1,5);
  vy_vec(ii) = dd(1,6);
  vz_vec(ii) = dd(1,7);


  XX = dd(:, 2);
  YY = dd(:, 3);
  ZZ = dd(:, 4);

  hold off
  plot(XX, YY, 'ok', 'linewidth',2)
  xlim([-0.0005, 0.0005])
  ylim([-0.0005, 0.0005])
  
  pause(0.001)
end 


figure
Nt = size(XX,2)
for ii = 1:size(XX,2)

  hold off
  plot(XX(:,1:ii)', YY(:, 1:ii)', 'k', 'linewidth', 2)
  hold on
  plot(XX(:, ii), YY(:, ii), 'ok')

%   xlim([min(min(XX)), max(max(XX))])
%   ylim([min(min(YY)), max(max(YY))])

  xlim([-0.0005, 0.0005])
  ylim([-0.0005, 0.0005])

  pause(0.001)
end

figure 
plot(t_vec, vy_vec, '-+b', 'linewidth', 2)

figure
plot(vy_vec, vz_vec, 'or', 'linewidth', 2)

figure
hold off
for ii = 1:1:numel(files_list)

  hold off
  plot(vx_vec(1:ii), vy_vec(1:ii), 'r', 'linewidth', 2)
  hold on
  plot(vx_vec(ii), vy_vec(ii), 'or', 'linewidth', 2)
  pause(0.001)
end


