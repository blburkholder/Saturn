function [q_MHD_strong_PS, q_MHD_strong] = get_heating_coefficient_MHD_strong_PS(power_spectrum, f, density, k_perp)
    vacum_permeability = 4*pi()*1e-7; %N/A^2
  
    del_B_3 = (power_spectrum.*f).^(3/2);
   
    q_MHD_strong_PS = 1/sqrt(density*vacum_permeability^3)*del_B_3.*k_perp; 
    q_MHD_strong = mean(q_MHD_strong_PS);  
end