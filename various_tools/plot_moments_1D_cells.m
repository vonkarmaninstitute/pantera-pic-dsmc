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

  qx(:, ii) = dd(:, 15);
  qy(:, ii) = dd(:, 16);
  qz(:, ii) = dd(:, 17);

end 


figure
for ii = 1:size(nn, 2)
  subplot(2,3,1)
  hold off
  plot(x_cc, nn(:,ii), '+k', 'linewidth', 2)
  hold on
  title(num2str(ii))
  xlabel('x [m]')
  ylabel('Npart [-]')

  subplot(2,3,2)
  hold off
  plot(x_cc, Ux(:,ii), '+k', 'linewidth', 2)
  hold on
  title(num2str(ii))
  xlabel('x [m]')
  ylabel('Ux [m/s]')


  subplot(2,3,3)
  hold off
  plot(x_cc, Uy(:,ii), '+k', 'linewidth', 2)
  hold on
  title(num2str(ii))
  xlabel('x [m]')
  ylabel('Uy [m/s]')

  subplot(2,3,4)
  hold off
  plot(x_cc, qx(:,ii), '+k', 'linewidth', 2)
  hold on
  title(num2str(ii))
  xlabel('x [m]')
  ylabel('qx')


  subplot(2,3,5)
  hold off
  plot(x_cc, qy(:,ii), '+k', 'linewidth', 2)
  hold on
  title(num2str(ii))
  xlabel('x [m]')
  ylabel('qy')


  subplot(2,3,6)
  hold off
  plot(x_cc, qz(:,ii), '+k', 'linewidth', 2)
  hold on
  title(num2str(ii))
  xlabel('x [m]')
  ylabel('qz')


  pause(0.1)
end

figure
for ii = 1:size(nn, 2)
  hold off
  plot(x_cc, nn(:,ii), '+k', 'linewidth', 2)
  hold on
  plot(-x_cc, nn(:,ii), '--or', 'linewidth', 2)

  title(num2str(ii))
  xlabel('x [m]')
%   ylim([0, 400])
  pause(0.1)
end
