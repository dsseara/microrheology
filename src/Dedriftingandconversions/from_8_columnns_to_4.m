function[poss] = from_8_columnns_to_4(res)

% Very simple program, takes in the 8 columns out of trackmem and makes 
% them into the only 4 columns we really need, x- and y-position, frame
% number and bead ID. 
%
% INPUT : Matrix after trackmem
% OUTPUT : Same matrix with information from feature finding taken out.

poss=res(:,1:2);
poss(:,3)=res(:,6);
poss(:,4)=res(:,8);