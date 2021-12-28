close all
clear
clc

page_screen_output(0);

% Plots solution in time

files_list = dir('../dumps/sol_*');

figure
for ii = 1:1:numel(files_list)

  % Load file
  dd = load(['../dumps/',files_list(ii).name]);
  fprintf('Data from: %s\n', files_list(ii).name);

  % Extract data
  t_now = dd(1,1);  % first column is time
  xx = dd(:,2);     % position is in the second column
  UU = dd(:,3:end); % solution (primitive variables) is written in all other columns

  % Print stuff
  subplot(3,1,1)
  plot(xx, UU(:,1), 'linewidth', 2);
  hold off
  grid on
  title(['Time: ', num2str(t_now), ' s'])
  xlabel('Position [m]')
  ylabel('prim(1)')

  subplot(3,1,2)
  plot(xx, UU(:,2), 'linewidth', 2);
  hold off
  grid on
  xlabel('Position [m]')
  ylabel('prim(2)')

  subplot(3,1,3)
  plot(xx, UU(:,3), 'linewidth', 2);
  hold off
  grid on
  xlabel('Position [m]')
  ylabel('prim(3)')

  pause(0.002)

end


