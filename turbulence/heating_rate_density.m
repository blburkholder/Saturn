%ONLY FOR CASSINI DATA INSIDE MAGSPHERE
%BUT SIMULATION DATA WILL HAVE MASS = 1
function [q_MHD_strong,d_perp] = heating_rate_density(bx,by,bz,ux,uz,rhon,delta_t,ion_kg,H)

    ion_molar_mass = ion_kg*0.001;
    avagdro_number = 6.022140857e23;
    ion_mass = ion_molar_mass/avagdro_number;
    mu_0 = 4*pi()*1e-7; %N/A^2

    B_x = bx;
    B_y = by;
    B_z = bz;
    B_tot = sqrt(B_x.^2+B_y.^2+B_z.^2);
    B_total_mean = mean(B_tot);

    B_vector_mean = get_B_vector_mean(mean(bx), mean(by), mean(bz));

    [B_fluctuation_parallel, B_fluctuation_perp, B_std_parallel, B_std_perp] =...
        get_B_std_vector_components(B_vector_mean, B_x, B_y, B_z);

    [f, power_spectrum_perp] = get_power_spectrum_morlet(B_fluctuation_perp, B_std_perp, delta_t);
    [~, power_spectrum_par] = get_power_spectrum_morlet(B_fluctuation_parallel/B_total_mean, B_std_parallel, delta_t);
    %[f, power_spectrum_par] = get_power_spectrum_fft(B_fluctuation_parallel/B_total_mean, B_std_parallel, delta_t);
    rho = rhon*ion_mass;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    k_prependicular_MHD = get_k_prependicular(f, ux, uz, B_vector_mean);

    if length(f) > 0
        [q_MHD_strong_PS, q_MHD_strong] = get_heating_coefficient_MHD_strong_PS(...
            power_spectrum_perp,f,rho,k_prependicular_MHD);
    else
        q_MHD_strong = 0;
    end

    %Bf2 = mean(sum(B_fluctuation_parallel.^2))
    Bf2 = mean(k_prependicular_MHD.^2.*power_spectrum_par.^2);
    T = get_temperature(H);
    larmori = get_larmor_radius(T, sqrt(B_total_mean),.018);
    jays_int = evalintegral(k_prependicular_MHD*larmori);
    Va = sqrt(B_total_mean)./sqrt(mu_0*rho);
    kz = 2*pi/(60000000);
    d_perp = sqrt(pi/2)*Bf2*Va/(kz)*jays_int;
end
