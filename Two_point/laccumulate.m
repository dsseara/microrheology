
function [res]=laccumulate(list,binparams,dim)

% ;  accumulates total and total-squared results for the 'list' data
% ;  of get_corr, uses logarithmic binning and 'uniq' tricks to be fast.

%sumlist=sum(list,2)
%size(list)
%list
%; get the parameters
lrmin = binparams(1);
lrmax = binparams(2);
lrbinsize = binparams(3);
nbins = binparams(4);

%; declare the array we'll accumulate into
res = zeros((2*dim)+1,nbins);

% turn the r's into logarithmic bin numbers
r = floor( (log(list(1,:))-lrmin)/lrbinsize );
%sumr=sum(r)
nel = length( list(1,:) );

%; do the sorting and uniq'ing
[sr,s] = sort(r);
%sr = sr
slist = list(:,s);
[elsr,isr]=unique(sr);
u = [0,isr];
nu = length(u);
sr=sr+1;

for i=1:nu-1
    %i
    % Sadly, the following check is needed because of roundoff error.
    if sr(u(i)+1) > 0 && sr(u(i)+1) <= nbins
       %test = sr(u(i)+1)+1
       res(1:dim,sr(u(i)+1)) = sum( slist(2:dim+1,u(i)+1:u(i+1)) ,2);
       res(dim+1:2*dim,sr(u(i)+1)) = sum( slist(2:dim+1,u(i)+1:u(i+1)).^2 ,2);
       res(2*dim+1,sr(u(i)+1)) = u(i+1)-u(i);
    end
end

%res