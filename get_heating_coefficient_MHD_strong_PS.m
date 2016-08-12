function [q_MHD_strong_PS, q_MHD_strong] = get_heating_coefficient_MHD_strong_PS(B_total_std, B_total_mean, power_spectrum, f, density, l_prependicular, v_phi_rel, k_prependicular)
    vacum_permeability = 4*pi()*1e-7; %N/A^2
   
%     k_prependicular = 2*pi()*f/v_phi_rel;
    
%     k_prependicular_MHD./k_prependicular
    
    l_prependicular_scaling_constant = l_prependicular*k_prependicular(1)/(2*pi());
    
    l_prependicular_calc = l_prependicular_scaling_constant*2*pi()./k_prependicular;

    del_B_3 = (power_spectrum.*f).^(3/2);
    
    q_MHD_strong_PS = 1/sqrt(density*vacum_permeability^3)*del_B_3.*k_prependicular; 
    q_MHD_strong = mean(q_MHD_strong_PS);  
end