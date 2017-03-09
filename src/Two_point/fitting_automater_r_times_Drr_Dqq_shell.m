function[msdd_0 flor_0] = fitting_automater_r_times_Drr_Dqq_shell(data,tau,minradius,maxradius,displaying,beadradius,rDrrpath,savingplot)

% This function automates the fitting and plotting of r*Drr vs r. The fit is
% done with a logarithmic average, its result is A. Hence the two-point msd
% is calculated with : 
% msdD(tau) = 2 * A / a. 
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
% msdd - the two-point MSD calculated by the logarithmic average.


titletext = {['\tau = ' num2str(data(tau,1,2)) 's']};

r = (data(tau,:,1))';
Drr = (data(tau,:,3))';
Dqq = (data(tau,:,4))';

limpos = find(Drr>0&isfinite(Drr));
Drrpos = Drr(limpos);
rpos = r(limpos);

limqpos = find(Dqq>0&isfinite(Dqq));
Dqqpos = Dqq(limqpos);
rqpos = r(limqpos);

limneg = find(Drr<0);
Drrneg = Drr(limneg);
rneg = r(limneg);

Drrfil = Drrpos;
rfil = rpos;
Dqqfil=Dqqpos;
rqfil=rqpos;


limfil = find(rpos<maxradius & rpos>minradius);
Drrfilpos = Drrfil(limfil);
rfilpos = rfil(limfil);
limqfil = find(rqpos<maxradius & rqpos>minradius);
Dqqfilpos = Dqqfil(limqfil);
rqfilpos = rqfil(limqfil);

msd2p=data(tau,limpos, : );
msd2p=msd2p(:,limfil, : );
[msdd_0 flor_0]=msddE( msd2p,minradius,maxradius,0,beadradius,0);
msd2pq=data(tau,limqpos, : );
msd2pq=msd2pq(:,limqfil, : );
[msddq_0 florq_0]=msddE( msd2pq,minradius,maxradius,0,beadradius,0);

msdd_0(:,[3 8])=msddq_0(:,[3 8]);
A=msdd_0(2)*3*beadradius/2;
B=msdd_0(3)*3*beadradius/4;

% display error bars
rerr=r;
Drrerr=abs(Drr);
errDrr=sqrt(data(tau,:,5)./max(data(tau,:,7),1));
errrDrr=rerr.*errDrr';

figure
loglog(rpos,rpos.*Drrpos,'bo', 'LineWidth', 1)
hold on
loglog(rfilpos,rfilpos.*Drrfilpos,'ro', 'LineWidth', 1)
plot(rneg,-rneg.*Drrneg,'co')
line([1 1000],[A A],'color','k', 'Linewidth', 1)
plot(10 ,0.00001,'w.')

errorbar( rerr, rerr.*Drrerr, min(errrDrr-rerr.*Drrerr,-1e-10)+rerr.*Drrerr, errrDrr, 'k', 'LineStyle', 'none' )
errorbarlogx

xlabel('Separation distance r (\mum)','fontsize',26)
ylabel('r*D_{rr}(\mum^3)','fontsize',26)
set(gca,'fontsize',24, 'XTick',[1 10 100 1000],...
  'TickLength',[0.025 0.06],...
  'YTick',[1e-7 1e-6 1e-5 0.0001 0.001 0.01 0.1 1 10 100 1000 1e4 1e5 1e6 1e7] );
ymin=min(min(abs(data(:,:,1).*data(:,:,3))));
ymax=max(max(data(:,:,1).*data(:,:,3)));
ylim([min(0.00001, 10^floor(log10(ymin))) max(1, 10^ceil(log10(ymax)))])
xlim( [1 1000] );
title(titletext,'fontsize',26)

if savingplot == 1
    print('-dpng','-r150', [rDrrpath 'r_times_Drr_tau_' num2str(tau) '.png'])
    saveas(gcf,[rDrrpath 'r_times_Drr_tau_' num2str(tau) '.fig'])
    close all
end

%%%%%%%%%%%%%%%%%%%%5
% Dqq plot
limqneg = find(Dqq<0);
Dqqneg = Dqq(limqneg);
rqneg = r(limqneg);

% display error bars
rqerr=r;
Dqqerr=abs(Dqq);
errDqq=sqrt(data(tau,:,6)./max(data(tau,:,7),1));
errrDqq=rqerr.*errDqq';

figure
loglog(rqpos,rqpos.*Dqqpos,'bo', 'LineWidth', 1)
hold on
loglog(rqfilpos,rqfilpos.*Dqqfilpos,'ro', 'LineWidth', 1)
plot(rqneg,-rqneg.*Dqqneg,'co')
line([1 1000],[B B],'color','k', 'Linewidth', 1)
plot(10 ,0.00001,'w.')

errorbar( rqerr, rqerr.*Dqqerr, min(errrDqq-rqerr.*Dqqerr,-1e-10)+rqerr.*Dqqerr, errrDqq, 'k', 'LineStyle', 'none' )
errorbarlogx

xlabel('Separation distance r (\mum)','fontsize',26)
ylabel('r*D_{\theta\theta}(\mum^3)','fontsize',26)
set(gca,'fontsize',24, 'XTick',[1 10 100 1000],...
  'TickLength',[0.025 0.06],...
  'YTick',[1e-7 1e-6 1e-5 0.0001 0.001 0.01 0.1 1 10 100 1000 1e4 1e5 1e6 1e7] );
ylim([min(0.00001, 10^floor(log10(ymin))) max(1, 10^ceil(log10(ymax)))])
xlim( [1 1000] );
title(titletext,'fontsize',26)



if displaying == 1
disp([(data(tau,:,3))',(1:length(data(1,:,3)))'])
end


if savingplot == 1
    print('-dpng','-r150', [rDrrpath 'r_times_Dqq_tau_' num2str(tau) '.png'])
    saveas(gcf,[rDrrpath 'r_times_Dqq_tau_' num2str(tau) '.fig'])
    close all
end