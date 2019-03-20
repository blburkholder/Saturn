function [e1,e2,e3] = B_var(bxs,bys,bzs)
    bx_m = mean(bxs);
    by_m = mean(bys);
    bz_m = mean(bzs);
    bx_m2 = bx_m^2;
    bx_2m = mean(bxs.^2);
    by_m2 = by_m^2;
    by_2m = mean(bys.^2);
    bz_m2 = bz_m^2;
    bz_2m = mean(bzs.^2);
    bxz_m = mean(bxs.*bzs);
    bxy_m = mean(bxs.*bys);
    byz_m = mean(bys.*bzs);

    MB = [ bx_2m - bx_m2, bxy_m-bx_m*by_m, bxz_m-bx_m*bz_m;...
        bxy_m-bx_m*by_m, by_2m-by_m2, byz_m-by_m*bz_m;...
        bxz_m-bx_m*bz_m, byz_m-by_m*bz_m, bz_2m-bz_m2];
    [~,e] = eig(MB);
    e1 = e(1,1);
    e2 = e(2,2);
    e3 = e(3,3);
end
