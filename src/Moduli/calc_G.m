% calling command: [ omega,Gs,Gp,Gpp, dd, dda ] = calc_G(tau,msd,a,dim,T,clip,width)

%        ***Mason-Weitz micro-rheology in the world of IDL!***
% 
%     INPUTS:
%        tau: in seconds
%        msd: in microns^2
%        a: radius in microns
%        dim: the dimensionality of the data 
%        T: in Kelvin
%        clip: a fraction of G(s) below which G'(w) and G"(w) are
%        meaningless (see NOTES below)
%        width: the width of the gaussian that is used for the polyfit (see Notes below)
%
%     NOTES:
%
%        It needs more than a 7-8 points per decade of time/frequency
%        and really hates noise and long-wavelength ripples in the data.
% 
%        G'(w) and G"(w) are clipped at 0.03x G(s) as they are almost
%        certainly meaningless below that point unless the data is
%        *extremely* clean.  Set 'clip' to less than 0.03 to see more.
%        See Tom Mason's paper: PRL *79*, 3284, (1997) for details.
%
%        set the width to something bigger for noisy data, aok if the data
%        is very weakly curved. recommended starting value: 0.7
% 
%     OUTPUTS:
%                omega : is the frequency (s or omega), sec^{-1}.
%                Gs : is G(s)  in Pascals
%                Gp : is G'(w) in Pascals
%                Gpp : is G"(w) in Pascals
%                dd : second log derivative of the MSD
%                dda : second log derivative of Gs
%        Remember: in MKS (Pascals), water = 0.001 Pascal*secs = 0.01 Poise

%     REVISION HISTORY:
%                Written by John C. Crocker, U. Penn:    5/'99
%                Made 'width' a keyword, JCC             7/'99
%                Change Gp,Gpp formulae, JCC             7/'99
% 
%        This code 'micrheo.pro' is copyright 1999 by John C. Crocker
%        and may be freely distributed and used if properly attributed.

function [ omega,Gs,Gp,Gpp,dd,dda ] = calc_G(tau,msd,a,dim,T,clip,width)

tau=(tau(:))' ;
msd=(msd(:))' ;

% set up the 'constants'
kB = 1.38065e-23;                        % MKS
am = a*1e-6;                             % convert microns to meters
dt = tau;
omega = 1./dt;
msdm = msd*1e-12;                        % convert msd to meters
C = dim*kB*T/(3*pi*am);                  % multiply by the dimensionality/3.
foo = (pi/2)-1;                          % a handy constant
% if not keyword_set(clip) $
%        then clip = 0.03                ; throw away 1.5 decades down

% use 2nd order local formula for G(s)-- good to 1% of G(s)
[m,d,dd] = logderive(dt,msdm,width);
Gs = C./((m.*gamma(1+d)).*(1 +(dd/2)));

% use 2nd order local formula for G'(w),G"(w)-- good to 2% of G(w)
[g,da,dda] = logderive(omega,Gs,width);
%Gp  = g.*(1./(1+dda)).*  cos( (pi/2)*da - foo*da.*dda );       %as in the paper
%Gpp = g.*(1./(1+dda)).*  sin( (pi/2)*da - foo*(1-da).*dda );   %as in the paper
Gp  = g.*(1./(1+dda)).* ( cos( (pi/2)*da ) - foo*da.*dda );      %as in the code
Gpp = g.*(1./(1+dda)).* ( sin( (pi/2)*da ) - foo*(1-da).*dda );    %as in the code

% clip off the suspicious (i.e. unreliable) data
w = find(Gp < Gs*clip); nw=size(w);
if nw > 0 , Gp(w)=0; end
w = find(Gpp < Gs*clip); nw=size(w);
if nw > 0 , Gpp(w)=0; end

if (max(abs(dd)) > 0.15) || (max(abs(dda)) > 0.15)   %original value: 0.15
    warning ('Warning, high curvature in data, moduli may be unreliable!')
end