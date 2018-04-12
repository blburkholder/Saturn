active_date = starting_active_time;

ion_molar_mass = 18.01528e-3; %kg/mol (water)
avagdro_number = 6.022140857e23;    
ion_mass = ion_molar_mass/avagdro_number;

qks = zeros(length(active_date),1);
qms = zeros(length(active_date),1);

i = 2;
[location_data] = get_location_data(); 
[magnetometer_data] = get_magnetometer_data(i);
[n, length_of_magnetometer_data] = size(magnetometer_data);
dates = 24*60*(datenum(magnetometer_data(1,:)...
    , magnetometer_data(2,:), magnetometer_data(3,:)...
    , magnetometer_data(4,:), magnetometer_data(5,:)...
    , floor(magnetometer_data(6,:))));

%for gg = 1:length(active_date)
for gg = 1:length(active_date)
    while active_date(gg) > dates(end)
        i = i +1;
        [magnetometer_data] = get_magnetometer_data(i);
        [n, length_of_magnetometer_data] = size(magnetometer_data);

        dates = 24*60*(datenum(magnetometer_data(1,:)...
            , magnetometer_data(2,:), magnetometer_data(3,:)...
            , magnetometer_data(4,:), magnetometer_data(5,:)...
            , floor(magnetometer_data(6,:))));
    end

    ind = find(dates == active_date(gg));
    if ~isempty(ind) && ind+360 < length(magnetometer_data)
        X = magnetometer_data(7,ind:ind+360);
        Y = magnetometer_data(8,ind:ind+360);
        Z = magnetometer_data(9,ind:ind+360);

        location_initial = location_data([8,11], location_data(1,:) == magnetometer_data(1,ind) ...
            & location_data(2,:) == magnetometer_data(2,ind) & location_data(3,:) == magnetometer_data(3,ind) ...
            & location_data(4,:) == magnetometer_data(4,ind));

%     if ~isempty(ind) && ind+180 < length(magnetometer_data) && ind-180 > 0
%         X = magnetometer_data(7,ind-180:ind+180);
%         Y = magnetometer_data(8,ind-180:ind+180);
%         Z = magnetometer_data(9,ind-180:ind+180);
% 
%         location_initial = location_data([8,11], location_data(1,:) == magnetometer_data(1,ind-180) ...
%             & location_data(2,:) == magnetometer_data(2,ind-180) & location_data(3,:) == magnetometer_data(3,ind-180) ...
%             & location_data(4,:) == magnetometer_data(4,ind-180));

        B_r = X';
        B_theta = Y';
        B_phi = Z';
        B_tot = sqrt(X'.^2+Y'.^2+Z'.^2);

        N = length(B_r);
        time_step = 1; %1s
        lag = floor(N/10);

        B_r = B_r*1e-9;
        B_theta = B_theta*1e-9;
        B_phi = B_phi*1e-9;
        B_tot = B_tot*1e-9;   

        B_r_std = std(B_r);
        B_r_rms = rms(B_r);
        B_r_mean = mean(B_r);
        B_r_min = min(B_r);
        B_r_max = max(B_r);

        B_theta_std = std(B_theta);
        B_theta_rms = rms(B_theta);
        B_theta_mean = mean(B_theta);
        B_theta_min = min(B_theta);
        B_theta_max = max(B_theta);

        B_phi_std = std(B_phi);
        B_phi_rms = rms(B_phi);
        B_phi_mean = mean(B_phi);
        B_phi_min = min(B_phi);
        B_phi_max = max(B_phi);

        B_total_rms = rms(B_tot);
        B_total_mean = mean(B_tot);
        B_total_min = min(B_tot);
        B_total_max = max(B_tot);

        B_vector_mean = get_B_vector_mean(B_r_mean, B_theta_mean, B_phi_mean);
        [B_fluctuation_parallel, B_fluctuation_perp, B_std_parallel, B_std_perp] = get_B_std_vector_components(B_vector_mean, B_r, B_theta, B_phi);
        tau = 50; %s

        B_total_std = B_std_perp;
        delta_t = 1;

        [f, power_spectrum_morlet] = get_power_spectrum_morlet(B_fluctuation_perp, B_std_perp, delta_t);
        %[f, power_spectrum_fft] = get_power_spectrum_fft(B_fluctuation_perp, B_std_perp, delta_t);

        [gyration_frequency] = get_gyration_frequency(B_tot, 1);
        H = get_scale_height(mean(location_initial(2)));

        [v_phi_rel, v_r_rel] = get_v_rel(location_initial(2));
        number_density = get_density(location_initial(1), H, location_initial(2));
        density = number_density*ion_mass;
        T = get_temperature(H);

        larmor_radius = get_larmor_radius(T, B_total_mean, 1);

        l_prependicular = get_magnetic_field_prependicular_length(tau, v_phi_rel, v_r_rel, 1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        q_KAW_delB = get_heating_coefficient_KAW(B_total_std, larmor_radius, density, l_prependicular);
        q_MHD_strong_delB = get_heating_coefficient_MHD_strong(B_total_std, B_total_mean, density, l_prependicular);

        [coeff_MHD, coeff_KAW, f_MHD, f_KAW, power_spectrum_MHD, power_spectrum_KAW] = get_power_spectrum_slopes(f, power_spectrum_morlet, gyration_frequency);      
        %[coeff_MHD, coeff_KAW, f_MHD, f_KAW, power_spectrum_MHD, power_spectrum_KAW] = get_power_spectrum_slopes(f, power_spectrum_fft, gyration_frequency);       

        k_prependicular_MHD = get_k_prependicular(f_MHD, v_phi_rel, v_r_rel, B_vector_mean,1);
        k_prependicular_KAW = get_k_prependicular(f_KAW, v_phi_rel, v_r_rel, B_vector_mean,1);
        %Khi = get_Khi(B_total_std, B_total_mean, l_parallel, l_prependicular);

        if length(f_MHD) > 0
            [q_MHD_strong_PS, q_MHD_strong] = get_heating_coefficient_MHD_strong_PS(B_total_std, B_total_mean, power_spectrum_MHD, f_MHD, density, l_prependicular, v_phi_rel, k_prependicular_MHD);
        else
            q_MHD_strong = 0;
        end

        if  length(f_KAW) > 0
            [q_KAW_PS, q_KAW] = get_heating_coefficient_KAW_PS(B_total_std, power_spectrum_KAW, f_KAW, larmor_radius, density, k_prependicular_KAW);
        else
            q_KAW = 0;
        end  

        qks(gg) = q_KAW;
        qms(gg) = q_MHD_strong;
    end

%     figure
%     plot(X)
%     hold on
%     plot(Y)
%     plot(Z)
%     pause
end

xx = 1:length(active_date);
xxk = xx(qks>0);
qksm = qks(qks > 0);
xxm = xx(qms>0);
qmsm = qms(qms > 0);

figure
scatter(xxk,qksm,'b')
hold on
scatter(xxm,qmsm,'r')
set(gca,'yscale','log')

figure
hist(log10(qksm))
