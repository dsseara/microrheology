function [rg] = calc_rg_k(poscork,dim)

% This program calculates the radius of gyration for one section of a bead
% trajectory. It is mainly used as a subfunction, so trying to use it
% directly might not be a good idea.

    N=size(poscork,1);
    rgpart=0;
    for n=1:N
        for m=n:N
            if dim==3
                rgpart = rgpart + ( poscork(m,1) - poscork(n,1) )^2 + ( poscork(m,2) - poscork(n,2) )^2 + ( poscork(m,3) - poscork(n,3) )^2 ;
            elseif dim==2
                rgpart = rgpart + ( poscork(m,1) - poscork(n,1) )^2 + ( poscork(m,2) - poscork(n,2) )^2;
            end
        end
    end

    steps=max(poscork(:,3))-min(poscork(:,3));
    rg=sqrt(rgpart/steps^2);
    clear poscork;

    % There are 3 ways to define steps.
    % 1 - steps is the maximum frame number minus the minimum frame number 
    % steps = max(poscork(:,3))-min(poscork(:,3));
    % 2 - steps is the maximum number of frames (a bit sloppy)
    % steps = max(poscork(:,3));
    % 3 - steps is the trajectory time divided by the average time step
    % steps = max(bead(:,4))/av_time