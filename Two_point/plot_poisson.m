function plot_poisson( pois, min_w )%, x2, ymatrix2)
% Create figure
figure1 = figure('FileName','C:\Vincent\backup microscope cigale\poisson_.fig');

% Create axes
axes1 = axes(...
  'FontSize',24,...
  'TickLength',[0.025 0.06],...
  'XTick',[0.001 0.01 0.1 1 10 100],...
  'YTick',[-1.4 -1.2 -1 -.8 -.6 -.4 -.2 0 .2 .4 .6 .8 1],...
  'XScale','log',...
  'Parent',figure1);
axis([0.001 100 -1.5 1]);
box('on');
hold on
 
plot( 1./pois(1:min_w,1), pois(1:min_w,3), 'x', 'MarkerSize', 8 );
plot( [1e-3 100], [0.5 0.5], 'k--', [1e-3 100], [-1 -1], 'k--' )

disp( ['mean: ' num2str(mean(pois(1:min_w,3)))] )

% Create xlabel
xlabel('\omega (rad/s)','FontSize',24);
 
% Create ylabel
ylabel('Poisson ratio','FontSize',24);
 
