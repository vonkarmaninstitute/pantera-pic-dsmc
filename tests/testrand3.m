close all
clear
clc

dd = load('procnumbers2');

% Number of processors
NP = max(dd(:,1))+1;

% Initialize cell array
for ii = 1:NP
  data{ii} = [];
end

% Unpack data
for ii = 1:size(dd,1)

  pid = dd(ii,1) + 1;

  data{pid} = [data{pid}; dd(ii,2)];
end


colors = 'rgbmcykrgbmcykrgbmcykrgbmcykrgbmcykrgbmcykrgbmcykrgbmcykrgbmcyk'
figure
hold on
for ii = 1:NP
  plot(data{ii}(1:10), '-o', 'color', colors(ii))
end

