function [q_MHD] = no_vitaliy_bullshit_heating_rate_density(bx,by,bz,R,lat,delta_t)
    
    ion_molar_mass = 18;
    avagdro_number = 6.022140857e23;    
    ion_mass = ion_molar_mass/avagdro_number; %mass per an individual oxygen atom
    mu_0 = 4*pi()*1e-7; %N/A^2

    B_vector_mean = [mean(bx),mean(by),mean(bz)];

    length_B = length(bx);    
    B_vector = [bx, by, bz].';
    
    B_ave(1,1:length_B) = B_vector_mean(1);
    B_ave(2,1:length_B) = B_vector_mean(2);
    B_ave(3,1:length_B) = B_vector_mean(3);
    
    unit_B_ave = [(B_ave(1,:)./sqrt(B_ave(1,:).^2 + B_ave(2,:).^2 + B_ave(3,:).^2)).', ...
                 (B_ave(2,:)./sqrt(B_ave(1,:).^2 + B_ave(2,:).^2 + B_ave(3,:).^2)).', ...
                 (B_ave(3,:)./sqrt(B_ave(1,:).^2 + B_ave(2,:).^2 + B_ave(3,:).^2)).'].';
    
    B_pertr = B_vector - B_ave;
    
    B_parallel_pertr(1,:) = dot(B_pertr, unit_B_ave).*unit_B_ave(1,:);
    B_parallel_pertr(2,:) = dot(B_pertr, unit_B_ave).*unit_B_ave(2,:);
    B_parallel_pertr(3,:) = dot(B_pertr, unit_B_ave).*unit_B_ave(3,:);
    
    B_fluct_perp = B_pertr - B_parallel_pertr;
    
    NFFT = 2^(nextpow2(length(B_fluct_perp(1,:))));
    N = ceil((NFFT+1)/2);

    f_sample = 1/delta_t;
    f_perp = (0:N-1)*f_sample/NFFT;
    f_perp = f_perp(2:length(f_perp));

    power_spectrum_tot_x = sum(cwt(B_fluct_perp(1,:), 1./(1.03*f_perp), 'morl').^2, 2).'; 
    power_spectrum_tot_y = sum(cwt(B_fluct_perp(2,:), 1./(1.03*f_perp), 'morl').^2, 2).'; 
    power_spectrum_tot_z = sum(cwt(B_fluct_perp(3,:), 1./(1.03*f_perp), 'morl').^2, 2).'; 

    power_spectrum_tot = sqrt(power_spectrum_tot_x.^2 + power_spectrum_tot_y.^2 + power_spectrum_tot_z.^2);
    power_spectrum_perp = 2/N*power_spectrum_tot;
  
    del_B3_perp = (power_spectrum_perp.*f_perp).^(3/2);

    %empirical profiles from thomsen 2010
    H =sqrt(R^2/(3*8.7));
    uz = (7.2*R-23)*1e3; %this is v_phi
    ux = (-0.74*R+15)*1e3; %v_r
    n_0 = 51880e6*(1./R)^4.1;
    number_density = n_0*exp(-R^2*(1-cos(pi()/180*lat)^6)./(3*H^2));  
   
    v_rel = [ux, 0, uz];
    v_rel_mag = sqrt(v_rel(1)^2 + v_rel(2)^2 + v_rel(3)^2);
    B_vector_mean_mag = sqrt(B_vector_mean(1)^2 + B_vector_mean(2)^2 + B_vector_mean(3)^2);
    angle_VB = acos(dot(v_rel, B_vector_mean)/(v_rel_mag*B_vector_mean_mag)); 
    k_perp = 2*pi()*f_perp/(v_rel_mag*sin(angle_VB));

    q_MHD_strong_PS = 1/sqrt(ion_mass*number_density*mu_0^3)*del_B3_perp.*k_perp;

    %this is the integration step
    q_MHD = mean(q_MHD_strong_PS);  

end
      
            


              