close all
clear
clc

page_screen_output(0);

files_list = dir('../dumps/moments_*');
% files_list = dir('/media/starlight/Maxtor/PANTERA_data/proc_00000_time_*');
% files_list = dir('/media/starlight/Maxtor/PANTERA_data_test_Boris/proc_00000_time_*');

for ii = 1:numel(files_list)

  % Load file
  dd = load(['../dumps/',files_list(ii).name]);
  % dd = load(['/media/starlight/Maxtor/PANTERA_data_test_Boris/',files_list(ii).name]);
  fprintf('Data from: %s\n', files_list(ii).name)

  x_cc  = dd(:, 2);

  nn(:, ii) = dd(:, 5);

  Ux(:, ii) = dd(:, 6);
  Uy(:, ii) = dd(:, 7);
  Uz(:, ii) = dd(:, 8);

end 


figure
for ii = 1:size(nn, 2)
  hold off
  plot(x_cc, nn(:,ii), 'k', 'linewidth', 2)
  title(num2str(ii))
  xlabel('x [m]')
  ylim([0, 400])
  pause(0.1)
end
