% This function is based on matlab's polyfit function:

%POLYFIT Fit polynomial to data.
%   P = POLYFIT(X,Y,N) finds the coefficients of a polynomial P(X) of
%   degree N that fits the data Y best in a least-squares sense. P is a
%   row vector of length N+1 containing the polynomial coefficients in
%   descending powers, P(1)*X^N + P(2)*X^(N-1) +...+ P(N)*X + P(N+1).

% It was modified such that it has an aditional input w, which assignes a
% different weight to every data point. w needs to be of the same size as x
% and y. The fitting is now done giving each point the corresponding
% weight.

function [ p ] = polyfitw(x,y,w,n)

% make sure that x,y and w are of the same size
if ~isequal(size(x),size(y))
    error('MATLAB:polyfit:XYSizeMismatch',...
          'X and Y vectors must be the same size.')
end

if ~isequal(size(x),size(w))
    error('MATLAB:polyfit:XYSizeMismatch',...
          'X and W vectors must be the same size.')
end

% make sure that x,y and w are all column vectors.
x = x(:);
y = y(:);
w = w(:);

%constract w as a matrix
W=zeros(length(x),length(x));
for i=1:length(x)
    W(i,i)=w(i);
end
%W=W;

% Construct Vandermonde matrix
V(:,n+1) = ones(length(x),1,class(x));
for j = n:-1:1
   V(:,j) = x.*V(:,j+1);
end

%V=V;
coorm=(V'*W)*(V);
%croom=(w*V')*(V*w')
b=(V'*W) * (y);
coorm^-1;
% Solve least squares problem
%p = ((V'*W)*(W*V))^-1 * (V'*W) * (W*y);
p= coorm^-1 * b;

p = p.';          % Polynomial coefficients are row vectors by convention.