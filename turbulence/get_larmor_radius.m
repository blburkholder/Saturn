function larmor_radius=get_larmor_radius(T, B_mean,ion_kg)

    avagdro_number = 6.022140857e23; 
    k_B = 1.3806488e-23; %J/K
    Z = 1;
    e = 1.60217657e-19; %coulombs
    
    larmor_radius = sqrt(ion_kg/avagdro_number*k_B*T)/(Z*e*B_mean);

end
