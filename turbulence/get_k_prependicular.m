function [k_prependicular] = get_k_prependicular(f, v_phi_rel, v_r_rel, B_vector_mean)
    v_rel = [v_r_rel, 0, v_phi_rel];
    
    v_rel_mag = sqrt(v_rel(1)^2 + v_rel(2)^2 + v_rel(3)^2);
    
    B_vector_mean_mag = sqrt(B_vector_mean(1)^2 + B_vector_mean(2)^2 + B_vector_mean(3)^2);
    
    angle_VB = acos(dot(v_rel, B_vector_mean)/(v_rel_mag*B_vector_mean_mag));
    
    k_prependicular = 2*pi()*f/(v_rel_mag*sin(angle_VB));
end