function [msd1 flor] = msdd(msd2p,rmin,rmax,dtmax,a,width,lfit2)
%
%       A little function for making 'distinct' msd's out of D_par,D_perp
%
% if no a agument, then a = 1.0;
% if no width argument, width = 0.7;

% this function currently assumes array argument is array not filename

% get the parameters
dt = msd2p(:,1,2);
if dtmax ~= 0 %% in this instantiation of code, argument of dtmax=0 means place no upper limit on dt, include all taus
    dt = dt(find(dt < dtmax));
end
ndt = length(dt);
dr = msd2p(1,:,1);
w = find((dr > rmin) & (dr < rmax));
nr=length(w);
if nr==0 % Flow control
    error('No dr data in stated [rmin,rmax] interval!') % Error message display
end
dr = dr(w);

% calculate the msd_d's
msd2w = msd2p(:,w,:);
msd1 = zeros(ndt,6,'double');
msd1(:,1) = dt;
flor = zeros(ndt,2,'double');
for i=1:ndt
    par = msd2w(i,:,3);     % get the longitudinal mode
    perp = msd2w(i,:,4);    % get the transverse mode
    r = msd2w(i,:,1);       % get the radii
    rbar = sum(r)/nr;
    
    parr = par.*r;
    perpr = perp.*r;
    
    if lfit2 == 0  
        msd1(i,2) = (2/(3*a))*sum(parr)/length(parr);
        msd1(i,3) = (4/(3*a))*sum(perpr)/length(perpr);
    else
        fit1 = polyfit(r,parr,1);
        fit2 = polyfit(r,perpr,1);
        msd1(i,2) = fit1(2)*(2/(3*a));    % constant for longitudinal mode
        msd1(i,3) = fit2(2)*(4/(3*a));    % constant for transverse mode
        flor(i,1) = fit1(1)*(2/(3*a))*rbar;
        flor(i,2) = fit2(1)*(4/(3*a))*rbar;
    end
    
    % try a cheesy estimate of the errors, usually no good.
    err = sqrt(msd2w(i,:,5)./msd2w(i,:,7));
    msd1(i,6) = (2/3)*sqrt(sum((r.*err).^2)/length(par));
end

% make up smoother versions of msd's
% msd1(:,4) = lfit(msd1(:,1),msd1(:,2),width);
% msd1(:,5) = lfit(msd1(:,1),msd1(:,3),width);
% NOTE: done outside of msdd by call_3_fitting_methods, since that calls msdd with only one tau, so get lots of
% warnings about the fitting of only one point.
