function[] = number_of_points_per_r(data,tau,pathout,parttitle)

figure
semilogx(data(tau,:,1),data(tau,:,7),'.','color',[0 0 0],'markersize',24)

topone = max(data(1,:,7));

% ylim([0 topone*1.2])

titletext = {['Number of points per r, \tau = ' num2str(data(tau,1,2)) 's'];parttitle};

title(titletext,'fontsize',24)
xlabel('Separation distance r (\mum)','fontsize',24)
ylabel('Number of points','fontsize',24)
set( gca, 'FontSize', 24 )


saveas(gcf,[pathout 'Number_of_points_vs_radius\n_vs_r_tau_' num2str(tau) '.tif'])
saveas(gcf,[pathout 'Number_of_points_vs_radius\n_vs_r_tau_' num2str(tau) '.fig'])
close all