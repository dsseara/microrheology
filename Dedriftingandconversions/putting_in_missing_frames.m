function[fullres] = putting_in_missing_frames(res)

% This program is used to insert 'fake' frames to replace those that are
% missing after running trackmem; that is to say the frames that were used
% as memory. The program makes a weighted average of the frame before and
% the frame after the missing frame(s). Hence, the feature is assumed to
% have followed a straight line between those points. The program is made
% to accept values of memory that can exceed 1. 
%
% INPUT : The matrix of displacements you get from trackmem.
% OUTPUT : The matrix now filled with these 'fake' frames. 

dframe = circshift(res(:,3),-1) - res(:,3) - 1;
temp = find(dframe(:,1)>0);
difones(:,1:2) = [temp,dframe(temp)];

for i = 1:length(difones(:,1))
    shifter = difones(i,2);
    res(difones(i,1)+shifter+1:end+shifter,:) = res(difones(i,1)+1:end,:);
    for k = 1:shifter
        res(difones(i)+k,:) = ((shifter+1-k)*(res(difones(i),:))+(k*res(difones(i)+shifter+1,:)))/(shifter+1);
    end
    difones(:,1) = difones(:,1) + shifter;
end
fullres = res;