%%%% define global variables
function [j1] = eval_integral(Bfield,T_i,T_e,krhovec,Bf2,larmori)
global alpha lambda zeta lbess1 lbess2 n zetad

%%%% load parameters

    magnetic_field = Bfield; % in T
    t_ion = T_i; % in eV
    t_electron = T_e; % in eV
    %rhoi = 102*sqrt(t_ion)/magnetic_field; % in km
    rhoi = larmori;
    kpar = 2*pi/60000000; % in m
    %L_B = 730; % in km  %%%%%%%%%%%%%%% what is this?
    %rhoonL = rhoi/L_B; %%%%%%%%%%%%%%%% what is this?
    %rhoonL = 0.1;
    rhoonL = 0;
    kprho = kpar*rhoi;  % k_parallel rho_i
    tau = t_electron/t_ion; %T_e/T_i
    %beta = 4.03e-11*0.1*t_ion/(magnetic_field*1e-5)^2
    beta = 1;

%%%%%%%%% define result vectors

    int1 = [];
    int2 = [];
    int3 = [];
    int4 = [];
    int5 = [];
    int6 = [];

    d1vec = [];
    d2vec = [];
    d3vec = [];
    d4vec = [];
    d5vec = [];
    d6vec = [];

    zetavec = [];
    zetadvec = [];
    ind = 0;

%krhovec = [-5:0.01:5];
for krho = krhovec
    ind = ind+1;

    eta = tau*krho^2/(1 + (1 + tau)*krho^2);  %% phi = -eta psi
    lambda = sqrt(2)*krho;  %%% for argument of Bessel function

    %%% define omega/sqrt(2) kpar v_i from dispersion relation
    %%% omega = kpar V_A (1 + (1 + T_e/T_i) krho^2

    zeta0 = sqrt(1 + (1 + tau)*krho^2)/sqrt(beta);

    %%% define omega_d/sqrt(2) kpar v_i
    %%% omega_d = ky v_i rhoonL

    zetad = (krho/kprho)*rhoonL/sqrt(2); %%%%%%%% \omega_d/ ( k_|| v_ti sqrt(2) ) ??????
    
    zetas = -(2/beta)*zetad;
    
    zetae = tau*krho^2/sqrt(beta);
    
    c(1) = 1;
    c(2) = -zetas;
    c(3) = -zeta0^2;
    c(4) = -zetae^2*zetas;
    
    rts = roots(c);

    %%% include drifts
    %zetap = (-zetad*(-2/beta) + sqrt(zeta^2 + zetad^2*(4/beta^2)))/2;
    %zetam = (-zetad*(-2/beta) - sqrt(zeta^2 + zetad^2*(4/beta^2)))/2;

    %%% to check small krho approximation
    %zeta = 1/sqrt(beta);

    %%% to check for effect of drifts in dispersion relation
    %find(max(real(rts))
    %zeta = max(rts);
    zeta = zeta0;

    %%% to check for case with no magnetic drifts
    %zetad = 0;

    %%%%%%%%%%%%  compute the resonant particle integrals
    %%% I^0_{00}
    n=0;
    lbess1 = 0;
    lbess2 = 0;
    f1 = rombint('fun',0,5);

    %%% I^0_{40}
    n=4;
    lbess1 = 0;
    lbess2 = 0;
    f2 = rombint('fun',0,5);

    %%% I^2_{00}
    n=2;
    lbess1 = 0;
    lbess2 = 0;
    f3 = rombint('fun',0,5);

    %%% I^2_{11}
    n=2;
    lbess1 = 1;
    lbess2 = 1;
    f4 = rombint('fun',0,5);
    
    %%% I^1_{10}
    n=1;
    lbess1 = 1;
    lbess2 = 0;
    f5 = rombint('fun',0,5);
    
    %%% I^3_{10}
    n=3;
    lbess1 = 1;
    lbess2 = 0;
    f6 = rombint('fun',0,5);

    int1 = [int1 f1];
    int2 = [int2 f2];
    int3 = [int3 f3];
    int4 = [int4 f4];
    int5 = [int5 f5];
    int6 = [int6 f6];

    %%%% compute the diffusion coefficients

    alpha = zetad/zeta;  %%% ratio of omega_d/omega
    eta = (tau*2/beta)*krho*sqrt(Bf2(ind));
    d1vec  = [d1vec f1*eta*eta];
    d2vec = [d2vec f2*alpha^2*(1 + eta)^2];
    d3vec = [d3vec -2*f3*alpha*eta*(1 + eta)];

    %%%% for parallel magnetic field

    %ratio = sqrt(20); %%%% |E_TD/E_{||}|
    %sigma = ratio*eta;
    sigma = sqrt(Bf2(ind))*krho;

    d4vec = [d4vec 2*sigma*sigma*f4/krho^2]; 
    d5vec = [d5vec -2*sigma*eta*f5/krho];   
    d6vec = [d6vec 2*sigma*alpha*(1+eta)*f6/krho];

    %%%%% tabulate for dispersion relation

    zetavec = [zetavec zeta];
    zetadvec = [zetadvec zetad];
end


% figure
% subplot(221)
% plot(krhovec,int1,krhovec,int2,krhovec,int3,krhovec,int4,krhovec,int5,krhovec,int6)
% legend('I^0_{00}','I^4_{00}','I^2_{00}','I^2_{11}','I^1_{10}','I^3_{10}')
% subplot(222)
% plot(krhovec,d1vec,krhovec,d2vec,krhovec,d3vec,krhovec,d4vec,krhovec,d5vec,krhovec,d6vec)
% legend('d_1','d_2','d_3','d_4','d_5','d_6')
% dtot = d1vec + d2vec + d3vec + d4vec + d5vec + d6vec;
% subplot(223)
% plot(krhovec,dtot,krhovec,d4vec)
% legend('D_{tot}','D_{B_{||}}')
% 
% subplot(224)
% plot(krhovec,zetavec*sqrt(beta),krhovec,zetadvec)  %%% plot omega/k_par V_A

%print -dpdf -r300 nodrift
j1 = nansum(d1vec);
end





