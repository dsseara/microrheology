function[rdf2] = drift_loop(res,smoo)
   
    poss = from_8_columnns_to_4(res);  % This part makes the matrix of 4 
                                       % columns that is more convenient.
    
    
    
    drift=motion(poss, [1 0], 2);      % Finds the center of mass motion 
                                       % for each frame.
    npts=length(drift(:,1));
    t=drift(:,1);
    dr(:,1)=mot_eintegrate(drift(:,3));
    dr(:,2)=mot_eintegrate(drift(:,4));                                   

figure
plot( t/16, dr )%Time hardcoded (can be wrong) on the plot
xlabel('time (s)','fontsize',26)
ylabel('drift (blue x, green y) (\mum)','fontsize',26)
set(gca,'fontsize',24,...
'TickLength',[0.025 0.06]);
hold on
    for j = 1:2             % smooths the points smoo away from the end points, then fit a line
                            % to those endpoints to smooth them too.
     dr(:,j)=smooth_d(dr(:,j), smoo);
     smoo2=fix(smoo/2)+3;
     dr(:,j)=dr(:,j)';
     line1=polyfit(t(1:smoo2+1),dr(1:smoo2+1,j),1);%polyfit gives slope first and then intersection
     line2=polyfit(t(npts-smoo2+1:end), dr(npts-smoo2+1:end,j),1);%which is exactly the inverse of that from IDL
     dr(1:smoo2+1,j)=dr(smoo2+2,j)-line1(1)*t(smoo2+2)+line1(1)*t(1:smoo2+1);
     dr(npts-smoo2+1:end,j)=dr(npts-smoo2-1+1,j)-line2(1)*t(npts-smoo2-1+1)+line2(1)*t(npts-smoo2+1:end);
    end
    rdf=dedrift(poss, dr);
plot( t/16, dr )%Time hardcoded (can be wrong) on the plot    
    rdf2 = putting_in_missing_frames(rdf);
    
 

