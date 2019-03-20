function l_prependicular = get_magnetic_field_prependicular_length(tau, v_phi_rel, v_r_rel)
    v_rel = sqrt(v_phi_rel^2 + v_r_rel^2);
    l_prependicular = tau*v_rel;
end
