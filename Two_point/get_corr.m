%;   calculates 'brownian correlation' components for track data....
function [bres]=get_corr(t,lbinparams,dt,dim,dedrift)
% NOTE: t is trajectories, here (not time; time is in the ncol-1 column)


%; declare some arrays
bres = zeros((2*dim)+1,lbinparams(4));
ncol = length(t(1,:));
nid  = max(t(ncol));

%; get the 'r' ranges (squared).
rminsq = exp(2*lbinparams(1));
rmaxsq = exp(2*lbinparams(2));

% get the rows, we do all for dt le 3, and 'triple count' for bigger dts.
st = circshift(t,[-dt,0]);
%sumt=sum(t,2), sumst=sum(st,2)
if dt <= 3
    %     if keyword_set( erode ) then begin
    %        mask = fltarr((2*fix(erode))+1.)+1
    %        bad = dilate(reform(t(0,*) eq 0),mask)
    %        bad = bad or shift(bad,-dt)     ; oink!
    %        st2 = shift(t,0,erode)
    %        st3 = shift(t,0,-dt-erode)
    %        ;we require that the point not be in a dilated 'bad'
    %        ;mask (for gaps) nor within erode of an end.
    %        w = where( (st(ncol-2,*)-t(ncol-2,*)) eq dt and $
    %           (st(ncol-1,*)-t(ncol-1,*)) eq 0 and $
    %           (st2(ncol-1,*)-t(ncol-1,*)) eq 0 and $
    %           (st3(ncol-1,*)-t(ncol-1,*)) eq 0 and (not bad),nw)
    %     endif else begin
    
            % find all that: time t and time shifted by dt is equal to dt
            %                neither position (t and shifted by dt) is 0 (?)
            %                the bead id for t and shifted by t is the same
    w = find( (st(:,ncol-1)-t(:,ncol-1)) == dt &...
        (t(:,1) ~= 0) & (st(:,1) ~= 0) &...
        (st(:,ncol)-t(:,ncol)) == 0);
    nw=length(w);
    %endelse
else
    %     if keyword_set( erode ) then begin
    %        mask = fltarr((2*fix(erode))+1.)+1
    %        bad = dilate(reform(t(0,*) eq 0),mask)
    %        bad = bad or shift(bad,-dt)
    %        st2 = shift(t,0,erode)
    %        st3 = shift(t,0,-dt-erode)
    %        ;we require that the point not be in a dilated 'bad'
    %        ;mask (for gaps) nor within erode of an end.
    %        w = where( (st(ncol-2,*)-t(ncol-2,*)) eq dt and $
    %           (st(ncol-1,*)-t(ncol-1,*)) eq 0 and $
    %           (st2(ncol-1,*)-t(ncol-1,*)) eq 0 and $
    %           (st3(ncol-1,*)-t(ncol-1,*)) eq 0 and (not bad) and $
    %           ((t(ncol-2,*) mod dt eq 0) or $
    %            (t(ncol-2,*) mod dt eq fix(dt/3.)) or $
    %            (t(ncol-2,*) mod dt eq fix(2*dt/3.))) ,nw)
    %     endif else begin

            % find all that: time t and time shifted by dt is equal to dt
            %                neither position (t and shifted by dt) is 0 (?)
            %                the bead id for t and shifted by t is the same
            % with added 
    w = find( (st(:,ncol-1)-t(:,ncol-1)) == dt &...
        (st(:,ncol)-t(:,ncol)) == 0 &...
        (t(:,1) ~= 0) & (st(:,1) ~= 0) &...
        ( ( mod(t(:,ncol-1),dt) == 1 ) |...
        ( mod(t(:,ncol-1),dt) == floor(dt/3)+1 ) |...
        ( mod(t(:,ncol-1),dt) == floor(2*dt/3)+1 ) ) );
    nw=length(w);
    %sumw=sum(w)-nw
    %     endelse
end

%nw=nw, sumw=sum(w)

if nw > 1   % ; -------------------------- ERW

    %; make up the data         xs(1:dim) is the position at tof all beads
    % which are there at t + dt too
    % and xs(dim+1:2*dim) is the position of those beads at t + dt
    xs = zeros(nw,2*dim);
    xs(:,1:dim) = t(w,1:dim);
    xs(:,dim+1:2*dim) = st(w,1:dim);
    %sumxs=sum(sum(xs))
    %pxs=(xs(4,1:4))

    %; dedrift, if we're allowed to by the user.
    if dedrift==1
        dx = sum(xs(:,dim+1:2*dim)-xs(:,1:dim))./nw;
        for i=1:dim, xs(:,dim+i) = xs(:,dim+i) - dx(i); end;
    end

    %; get the time and id info
    ti = t(w,ncol-1);
    id = t(w,ncol);

    %; sort the data by time
    [tie,s] = sort(ti);
    
    %don't think this is important, just want to compare to IDL
%     if dt==153
%     uti=unique(tie) 
%     for i=1:length(uti),
%         fti=find(tie==uti(i));
%         if length(fti)>1, s(fti)=sort(s(fti),'descend'); end
%     end
%     s
%     end
    %if dt==153, ti, end
    xs = xs(s,:);               % akk sorted by time now
    ti = ti(s);
    id = id(s);
    %pxs=(xs(1,1:4))

    %; get the indices of the unique *times*
    [elu,u] = unique(ti);
    ntimes = length(u);
    u=u';
    u = [0,u];
    
    %; get the maximum number of beads in one frame pair
    su=circshift(u,[0,-1]);
    maxn = max( su(1:ntimes)-u(1:ntimes) );
    %maxn is a single number

    maxlist = maxn^2/2;
    nreport = 500;
    if maxn > 100, nreport = 250; end
    if maxn > 200, nreport = 50; end
    if maxn > 500, nreport = 10; end

    %; define a triangular raster scan list
    bamat=ones(maxn,maxn);
    for j = 1:maxn
        bamat(j,:) = bamat(j,:)*j;
    end
    bbmat=bamat';
    
    w = find(bbmat > bamat); nw0=length(w);
    %w12=w(1:2)
    if (nw0 > 0)
        bamat = bamat(w);
        bbmat = bbmat(w);
    end
                % these vectors will contain respectively 1 1 2 1 2 3 ...
                % and 2 3 3 4 4 4 .... so that they can be use as indices
                % to do pairs: 1&2, 1&3, 2&3, 1&4...
                
    %bamat15=bamat(1:5),bbmat15=bbmat(1:5)

    %; define the scratchpad (list)
    if maxlist > 2e4, buf = maxlist; else buf = 2e4; end   %; 2e4 element lists sort the fastest
    list = zeros(dim+1,buf);
    point = 0;

    %; loop over the times
    for i = 1:ntimes

        %; get the relevant data for the ith time. (all the xy positions at
        % ith time (first dim columns), and at ith time + dt (last dim columns)
        lxs = xs( u(i)+1:u(i+1),: );
        %if i==1, plxs=lxs(1,1), end
        ngood = length(lxs(:,1));

        if ngood > 1   %; we should check!

            %;     fast N^2 distance calculator, inspired by 'track'
            ntri = ngood*(ngood-1)/2;
            %if i==1, ntri=ntri, end
            amat = bamat(1:ntri);
            bmat = bbmat(1:ntri);
            %if i==1, pamat=amat(1:3), pbmat=bmat(1:3), end
            %if i==1, sum(sum(lxs)), end
            %size(lxs), max(bmat)
            rsq = sum( (lxs(amat,1:dim)-lxs(bmat,1:dim)).^2 ,2);   
                    % calculate the distance squared between 1&2, 1&3, 2&3,... 

            %if i==10, rsq=rsq, end

            w = find((rsq < rmaxsq) & (rsq > rminsq));
            nok=length(w);

            if nok > 0
                %if i==100, nok, end

                %;    calculate the sep. vectors
                r   = sqrt(rsq(w));         % distance between beads a & b at t=i (for binning)
                %if i==1,r(1:5),end
                amatw = amat(w);
                bmatw = bmat(w);
                rx  = lxs(amatw,1) - lxs(bmatw,1);  % distance in x between a & b at t=i
                ry  = lxs(amatw,2) - lxs(bmatw,2);
                xa1 = lxs(amatw,1);                 % position in x of a at t=i
                ya1 = lxs(amatw,2);
                xa2 = lxs(amatw,dim+1);             % position in x of a at t=i+dt
                ya2 = lxs(amatw,dim+2);
                xb1 = lxs(bmatw,1);
                yb1 = lxs(bmatw,2);
                xb2 = lxs(bmatw,dim+1);
                yb2 = lxs(bmatw,dim+2);
                dxa = xa2 - xa1;                    % displacement of a in x between t=i and t=i+dt
                dya = ya2 - ya1;
                dxb = xb2 - xb1;
                dyb = yb2 - yb1;

                if dim == 2

                    %; calculate the longitudinal part
                    rand('state', sum(100*clock));
                    ran = 1-2*(rand(1,nok) > 0.5);  %; randomize        %ran is a vector of 1 and -1
                    ran=ran';
                    nx  = rx.*ran./r;      %; unit vector r1     +- the x component of unit vector between a and b at t=i
                    ny  = ry.*ran./r;
                    ddl = ((dxa.*nx)+(dya.*ny)).*((dxb.*nx)+(dyb.*ny)); % component of displacement of a doted to unit vector above times disp of b
                                        % so "ran" in effect randomizes
                                        % whether r is from a to b or b to
                                        % a. I guess that helps with
                                        % numerical rounding errors.
                    %if i==1, ddl, end

                    %; calculate the transverse part
                    rand('state', sum(100*clock));
                    ran = 1-2*(rand(1,nok) > 0.5);  %; randomize        %ran is a vector of 1 and -1
                    ran=ran';
                    px  = ny.*ran;        %; ortho unit vector      switching the y->x and x->-y we get perp vector, randomize using ran
                    py  = -nx.*ran;
                    ddt = ((dxa.*px)+(dya.*py)).*((dxb.*px)+(dyb.*py)); % same dot product, but perp.
                    %if i==1 & dt==133, sizeddt=size(ddt), end

                    %; add to the list
                    list(1,point+1:point+nok) = r;
                    list(2,point+1:point+nok) = ddl;
                    list(3,point+1:point+nok) = ddt;
                    point = point+nok;
                    %if i==6, sum(list,2), end

                    %             else
                    %
                    %                 ; get a couple more sep vectors
                    %                 rz  = lxs(2,amatw) - lxs(2,bmatw)
                    %                 za1 = reform(lxs(2,amatw))
                    %                 za2 = reform(lxs(dim+2,amatw))
                    %                 zb1 = reform(lxs(2,bmatw))
                    %                 zb2 = reform(lxs(dim+2,bmatw))
                    %                 dza = za2 - za1
                    %                 dzb = zb2 - zb1
                    %
                    %                 ; calculate the longitudinal part
                    %                 ran = 1-2*(randomu(seed,nok) ge 0.5)  ; randomize
                    %                 nx  = rx*ran/r
                    %                 ny  = ry*ran/r
                    %                 nz  = rz*ran/r
                    %                 ddl = ((dxa*nx)+(dya*ny)+(dza*nz))* $
                    %                 ((dxb*nx)+(dyb*ny)+(dzb*nz))
                    %
                    %                 ; calculate the phi transverse part
                    %                 ran = 1-2*(randomu(seed,nok) ge 0.5)  ; randomize
                    %                 rxy = sqrt(rx^2 + ry^2)
                    %                 px1  = ry*ran/rxy
                    %                 py1  = -rx*ran/rxy
                    %                 ddp = ((dxa*px1)+(dya*py1))*((dxb*px1)+(dyb*py1))
                    %
                    %                 ; calculate the theta transverse part
                    %                 ran = 1-2*(randomu(seed,nok) ge 0.5)  ; randomize
                    %                 px2  = -(nz*py1)*ran
                    %                 py2  = (nz*px1)*ran
                    %                 pz2  = (nx*py1-ny*px1)*ran
                    %                 ddt = ((dxa*px2)+(dya*py2)+(dza*pz2))* $
                    %                 ((dxb*px2)+(dyb*py2)+(dzb*pz2))
                    %
                    %                 ; add to the list
                    %                 list(point+1:point+nok,0) = r
                    %                 list(point+1:point+nok,1) = ddl
                    %                 list(point+1:point+nok,2) = ddt
                    %                 list(point+1:point+nok,3) = ddp
                    %                 point = point+nok

                end
            end
        end
        %if i==1, alist=list(find(list~=0)), end
        
        %; do the running totals if the buffer gets full
        if point > (buf-maxlist)
            %sumbres_b=sum(bres,2)
            %point
            bres = bres + laccumulate(list(:,1:point),lbinparams,dim);
            %sum(bres(2*dim+1,:))
            point=0;
            %sumbres_a=sum(bres,2)
        end

        if mod(round(i),nreport)==0
            disp([num2str(round(i)),'   ',num2str(ntimes),'   ',num2str(sum(bres(2*dim+1,:))+point-1)]);
        end

    end

    %; finish the running totals
    if point > 0
        bres = bres + laccumulate(list(:,1:point),lbinparams,dim);
    end
%test1=sum(bres(2*dim+1,:)), test2=point
disp([num2str(round(i)),'   ',num2str(ntimes),'   ',num2str(sum(bres(2*dim+1,:)))]);
%; -------------------------- ERW

end
