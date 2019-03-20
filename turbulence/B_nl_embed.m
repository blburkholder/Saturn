function [tot,br1,br2,bt1,bt2,bp1,bp2] = B_nl_embed(B_r,B_theta,B_phi)    
    br1 = B_r(1:2:end);
    br2 = B_r(2:2:end);    
    bt1 = B_theta(1:2:end);
    bt2 = B_theta(2:2:end);
    bp1 = B_phi(1:2:end);
    bp2 = B_phi(2:2:end);

    ints_t = 0;
    ints_r = 0;
    ints_p = 0;

    for  q = 1:(length(bt1)-1)
        xt1 = bt1(q); xr1 = br1(q); xp1 = bp1(q);
        yt1 = bt2(q); yr1 = br2(q); yp1 = bp2(q);
        xt2 = bt1(q+1); xr2 = br1(q+1); xp2 = bp1(q+1);
        yt2 = bt2(q+1); yr2 = br2(q+1); yp2 = bp2(q+1);

        for p = q:(length(bt1)-1)
            if p~= q
            nxt1 = bt1(p); nxr1 = br1(p); nxp1 = bp1(p);
            nyt1 = bt2(p); nyr1 = br2(p); nyp1 = bp2(p);
            nxt2 = bt1(p+1); nxr2 = br1(p+1); nxp2 = bp1(p+1);
            nyt2 = bt2(p+1); nyr2 = br2(p+1); nyp2 = bp2(p+1);

                [xp,~] = polyxpoly([xt1,xt2],[yt1,yt2],[nxt1,nxt2],[nyt1,nyt2]);
                if ~isempty(xp)
                    ints_t = ints_t + 1;
                end
                [xp,~] = polyxpoly([xr1,xr2],[yr1,yr2],[nxr1,nxr2],[nyr1,nyr2]);
                if ~isempty(xp)
                    ints_r = ints_r + 1;
                end
                [xp,~] = polyxpoly([xp1,xp2],[yp1,yp2],[nxp1,nxp2],[nyp1,nyp2]);
                if ~isempty(xp)
                    ints_p = ints_p + 1;
                end
    %             plot(x1:(x2-x1)/10:x2,f)
    %             hold on
    %             scatter(xp,yp)
    %             plot(nx1:(nx2-nx1)/10:nx2,fn)
    %             hold off
    %             pause
            end
        end
    end
    s_t = abs(max(br1)-min(br1))+abs(max(br2)-min(br2));
    s_r = abs(max(bt1)-min(bt1))+abs(max(bt2)-min(bt2));
    s_p = abs(max(bp1)-min(bp1))+abs(max(bp2)-min(bp2));

    tot = (ints_t + ints_r + ints_p)*(s_t+s_r+s_p);