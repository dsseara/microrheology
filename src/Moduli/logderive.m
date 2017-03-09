function [f2,df,ddf] = logderive(x,f,width)

%       A special purpose routine for finding the first and
%       second logarithmic derivatives of slowly varying,
%       but noisy data.  It returns 'f2'-- a smoother version
%       of 'f' and the first and second log derivative of f2.

np = length(x);
df = zeros(1,np);
ddf = zeros(1,np);
f2 = zeros(1,np);
lx = log(x);
ly = log(f);

for i=1:np
    w = exp( -(lx-lx(i)).^2 / (2*width.^2) );         % a 'Gaussian'
    ww = find(w > 0.03);                           % truncate the gaussian, run faster
    %s(i)=length(ww);
    res = polyfitw(lx(ww),ly(ww),w(ww),2);
    f2(i) = exp(res(3) + res(2)*lx(i) + res(1)*(lx(i)^2));
    df(i) = res(2)+(2*res(1)*lx(i));
    ddf(i) = 2*res(1);
end

% figure,plot(lx,ly,'.')
% figure,plot(lx,df,'r.')
% figure,plot(lx,ddf,'g.')
% figure,plot(lx,s,'k')