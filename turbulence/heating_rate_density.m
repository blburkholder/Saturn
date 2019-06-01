%ONLY FOR CASSINI DATA INSIDE MAGSPHERE
%BUT SIMULATION DATA WILL HAVE MASS = 1
function [q_MHD_strong,dci,dbe] = heating_rate_density(bx,by,bz,ux,uz,rhon,delta_t,ion_kg,H)

    ion_molar_mass = ion_kg*0.001;
    avagdro_number = 6.022140857e23;
    ion_mass = ion_molar_mass/avagdro_number;
    mu_0 = 4*pi()*1e-7; %N/A^2
    kB = 1.3806488e-23; %J/K

    B_x = bx;
    B_y = by;
    B_z = bz;
    B_tot = sqrt(B_x.^2+B_y.^2+B_z.^2);
    B_total_mean = mean(B_tot);

    B_vector_mean = get_B_vector_mean(mean(bx), mean(by), mean(bz));

    [~, B_fluctuation_perp,~, B_std_perp] =...
        get_B_std_vector_components(B_vector_mean,B_x,B_y,B_z);
    [B_fluctuation_parallel, ~, B_std_parallel, ~] =...
        get_B_std_vector_components(B_vector_mean, B_x, B_y, B_z);

    [f_morl, power_spectrum_perp] = get_power_spectrum_morlet(B_fluctuation_perp, B_std_perp, delta_t);
    [f_fft, power_spectrum_par] = get_power_spectrum_fft(B_fluctuation_parallel/B_total_mean, B_std_parallel, delta_t);
    rho = rhon*ion_mass;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gyro_freak = get_gyration_frequency(B_tot,ion_molar_mass);
    k_prependicular_morl = get_k_prependicular(f_morl(f_morl > 2e-3 & f_morl < 2e-1 & f_morl < gyro_freak/2), ux, uz, B_vector_mean);
    k_prependicular_fft = get_k_prependicular(f_fft, ux, uz, B_vector_mean);

    [~, q_MHD_strong] = get_heating_coefficient_MHD_strong_PS(...
        power_spectrum_perp(f_morl > 2e-3 & f_morl < 2e-1 & f_morl < gyro_freak/2),f_morl(f_morl > 2e-3 & f_morl < 2e-1 & f_morl < gyro_freak/2),rho,k_prependicular_morl);

    Bf2 = f_fft.*power_spectrum_par;
    T_i = get_temperature(H);
    Ti = 500;
    Te = 200;
    larmori = get_larmor_radius(T_i, B_total_mean,0.018);
    Va = B_total_mean/sqrt(mu_0*rho);
    krho = k_prependicular_fft*larmori;
    %[j1,j4] = evalintegral(B_total_mean,Ti,Te,krho,Bf2,larmori);
    j1 = evalintegral(B_total_mean,Ti,Te,krho(krho <= 1),Bf2(krho < 1),larmori);
    kz = 2*pi/(60000000);
    dci = sqrt(pi/8)/kz*Va*j1;
%    dbe = sqrt(pi/8)/kz*Va*mean(Bf2)*exp(-1);
    dbe = sqrt(pi/8)/kz*Va*mean(krho(krho <= 1).^2.*Bf2(krho < 1))*exp(-1);
end
