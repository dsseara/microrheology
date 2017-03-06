function[msdd_1 flor_1] = fitting_automater_fit_to_line_shell(data,tau,minradius,maxradius,displaying,beadradius,linfitpath,savingplot)

% This function automates the fitting and plotting of Drr vs r. The fit is
% done by fitting r*Drr to a straight line and using the constant part as r*Drr: 
% msdD(tau) = 2 * cst /3a.
%
% INPUTS :
% data - this is the output of the twopoint function.
% tau - this is the value of tau to try (it is an index for the data
% matrix, and so should be an integer.)
% badones - this allows you to eliminate bad data points quickly, based on
% their index, the next parameter helps out for this
% displaying - a logical for displaying Drr and the associated index. The
% plot will only show positive values, hence you can associate the points
% you want to associate to their index. There is no need to eliminate
% negative values, these get taken care of in the code.
% beadradius - the radius (in microns) of the beads.
%
% OUTPUT:
% msdd - the two-point MSD calculated by the intercept.

titletext = {['\tau = ' num2str(data(tau,1,2)) 's']};

r = (data(tau,:,1))';
Drr = (data(tau,:,3))';

limpos = find(Drr>0&isfinite(Drr));
Drrpos = Drr(limpos);
rpos = r(limpos);

limneg = find(Drr<0);
Drrneg = Drr(limneg);
rneg = r(limneg);

Drrfil = Drrpos;
rfil = rpos;

limfil = find(rpos<maxradius & rpos>minradius);

Drrfilpos = Drrfil(limfil);
rfilpos = rfil(limfil);

msd2p=data(tau,limpos, : );
msd2p=msd2p(:,limfil, : );
[msdd_1 flor_1]=msddE( msd2p,minradius,maxradius,0,beadradius,1);

allx = 1:1000;

rbar = mean(rfilpos);
fity = (msdd_1(2) +flor_1(1)*allx/rbar)*3/2*beadradius;
fitpos=find(fity>0);

[msdd_0 flor_0]=msddE( msd2p,minradius,maxradius,0,beadradius,0);
A=msdd_0(2)*3*beadradius/2;


figure
loglog(rpos,Drrpos.*rpos,'bo')
hold on
line([1 1000],[A A],'color','k', 'Linewidth', 1, 'LineStyle', '--')
loglog(rfilpos,Drrfilpos.*rfilpos,'ro')
loglog(allx(fitpos),fity(fitpos),'k-', 'Linewidth', 1);
plot(rneg,-rneg.*Drrneg,'co')

plot(10 ,0.00001,'w.')
set( gca, 'XTick',[1 10 100 1000],...
  'TickLength',[0.025 0.06],...
  'YTick',[1e-7 1e-6 1e-5 0.0001 0.001 0.01 0.1 1 10 100 1000 1e4 1e5 1e6 1e7] );

ymin=min(min(abs(data(:,:,1).*data(:,:,3))));
ymax=max(max(data(:,:,1).*data(:,:,3)));
ylim([min(0.00001, 10^floor(log10(ymin))) max(1, 10^ceil(log10(ymax)))])
xlim( [1 1000] );

xlabel('Separation distance r (\mum)','fontsize',26)
ylabel('r*D_{rr}(\mum^3)','fontsize',26)
set(gca,'fontsize',24)
title(titletext,'fontsize',26)


if displaying == 1
disp([(data(tau,:,3))',(1:length(data(1,:,3)))'])
end


if savingplot == 1
    print('-dpng','-r150', [linfitpath '\rDrr_vs_r_fit_to_line_tau_' num2str(tau) '.png'])
    saveas(gcf,[linfitpath '\rDrr_vs_r_fit_to_line_tau_' num2str(tau) '.fig'])
    close all
end
