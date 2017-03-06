function f2 = lfit(x,f,width)
% based on lfit of John Crocker.
%       A special purpose routine for finding the first and
%       second logarithmic derivatives of slowly varying,
%       but noisy data.  It returns 'f2'-- a smoother version 
%       of 'f'. based on 'logderiv' in micrheo.
%   IDL version sets default width to 0.7.  Note:  this results in harsh smoothing!

% Smooths data which is logarithmic in x and y; will produce complex
% results if the y data has negative points and issue a warning (added
% by Vincent Pelletier on 07/1/7)

np = prod(size(x));
df = zeros(1,np);
ddf = zeros(1,np);
f2 = zeros(1,np);
lx = log(x);
ly = log(f);

fneg=find( f < 0 );
if length( fneg ) > 0
    warning( 'Negative value found, lfit will produce complex output.' );
end

for i=1:np
        w = exp( -((lx-lx(i)).^2) / (2.0*width^2) ); % a 'Gaussian'
        ww = find(w > 0.03);   % truncate the gaussian, run faster
        res = polyfit(lx(ww),ly(ww),2);  %!!!!!! no weighting yet
        % Matlab orders the polynomials in decreasing order, IDL increasing
        f2(i) = exp(res(3) + res(2)*lx(i) + res(1)*(lx(i)^2));
        df(i) = res(2)+(2.0*res(1)*lx(i));
        ddf(i) = 2.0*res(1);
end
