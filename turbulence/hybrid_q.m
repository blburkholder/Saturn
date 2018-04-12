x = ans.X;
bx = ans.BX;
bz = ans.BZ;
ux = ans.UX;
uz = ans.UZ;
rhon = ans.DENS;

%ion_molar_mass = 18.01528e-3; %kg/mol (water)
ion_molar_mass = 1e-3; %kg/mol (protons)
avagdro_number = 6.022140857e23;    
ion_mass = ion_molar_mass/avagdro_number;

B_x = bx;
B_y = zeros(size(bx));
B_z = bz;
B_tot = sqrt(B_x.^2+B_y.^2+B_z.^2);
  
B_x_std = std(B_x);
B_x_rms = rms(B_x);
B_x_mean = mean(B_x);
B_x_min = min(B_x);
B_x_max = max(B_x);

B_y_std = std(B_y);
B_y_rms = rms(B_y);
B_y_mean = mean(B_y);
B_y_min = min(B_y);
B_y_max = max(B_y);

B_z_std = std(B_z);
B_z_rms = rms(B_z);
B_z_mean = mean(B_z);
B_z_min = min(B_z);
B_z_max = max(B_z);

B_total_rms = rms(B_tot);
B_total_mean = mean(B_tot);
B_total_min = min(B_tot);
B_total_max = max(B_tot);

B_vector_mean = get_B_vector_mean(B_x_mean, B_y_mean, B_z_mean);
[B_fluctuation_parallel, B_fluctuation_perp, B_std_parallel, B_std_perp] = get_B_std_vector_components(B_vector_mean, B_x, B_y, B_z);

B_total_std = B_std_perp;
delta_t = 1;

[f, power_spectrum_morlet] = get_power_spectrum_morlet(B_fluctuation_perp, B_std_perp, delta_t);
%[f, power_spectrum_fft] = get_power_spectrum_fft(B_fluctuation_perp, B_std_perp, delta_t);

[gyration_frequency] = get_gyration_frequency(B_tot, 0);
T = 5e6;
tau = 50;
larmor_radius = get_larmor_radius(T, B_total_mean, 0);
l_prependicular = get_magnetic_field_prependicular_length(tau, mean(ux), mean(uz), 0);
rho = rhon*ion_mass;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[coeff_MHD, coeff_KAW, f_MHD, f_KAW, power_spectrum_MHD, power_spectrum_KAW] = get_power_spectrum_slopes(f, power_spectrum_morlet, gyration_frequency);      
%[coeff_MHD, coeff_KAW, f_MHD, f_KAW, power_spectrum_MHD, power_spectrum_KAW] = get_power_spectrum_slopes(f, power_spectrum_fft, gyration_frequency);       

k_prependicular_MHD = get_k_prependicular(f_MHD, mean(ux), mean(uz), B_vector_mean,1);
k_prependicular_KAW = get_k_prependicular(f_KAW, mean(ux), mean(uz), B_vector_mean,1);

if length(f_MHD) > 0
    [q_MHD_strong_PS, q_MHD_strong] = get_heating_coefficient_MHD_strong_PS(B_total_std, B_total_mean, power_spectrum_MHD, f_MHD, mean(rho), l_prependicular, ux, k_prependicular_MHD);
else
    q_MHD_strong = 0;
end

if  length(f_KAW) > 0
    [q_KAW_PS, q_KAW] = get_heating_coefficient_KAW_PS(B_total_std, power_spectrum_KAW, f_KAW, larmor_radius, mean(rho), k_prependicular_KAW);
else
    q_KAW = 0;
end 
 
q_MHD_strong
q_KAW
