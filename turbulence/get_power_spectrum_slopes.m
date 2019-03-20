function [coeff_MHD, coeff_KAW, f_MHD, f_KAW, power_spectrum_MHD, power_spectrum_KAW] = get_power_spectrum_slopes(f, power_spectrum, gyration_frequency,sim_or_cass)
    %sim_or_cass = 0 for simulation
    %sim_or_cass = 1 for cassini data   

%     figure
%     plot(log10(f),log10(power_spectrum))
%     hold on
%     a = [ones(length(power_spectrum),1),log10(f)']\...
%         log10(power_spectrum)';
%     plot([log10(f(1)),log10(f(end))],[a(1)+a(2)*log10(f(1)),a(1)+a(2)*log10(f(end))])
%     title(['slope = ',num2str(a(2))])
 
    if sim_or_cass
        %f_MHD = f(f < gyration_frequency/2 & f > 2E-3 & f < 2E-1);
        %power_spectrum_MHD = power_spectrum(f < gyration_frequency/2 & f > 2E-3 & f < 2E-1);
        f_MHD = f;
        power_spectrum_MHD = power_spectrum;

        f_KAW = f(f > gyration_frequency*2 & f < 2E-1 & f > 2E-3);
        power_spectrum_KAW = power_spectrum(f > gyration_frequency*2 & f < 2E-1 & f > 2E-3);  
    else
        f_MHD = f(5:end-20);
        power_spectrum_MHD = power_spectrum(5:end-20);

%         figure
%         plot(log10(f_MHD),log10(power_spectrum_MHD))
%         hold on
%         a = [ones(length(power_spectrum_MHD),1),log10(f_MHD)']\...
%             log10(power_spectrum_MHD)';
%         plot([log10(f_MHD(1)),log10(f_MHD(end))],[a(1)+a(2)*log10(f_MHD(1)),a(1)+a(2)*log10(f_MHD(end))])
%         title(['slope = ',num2str(a(2))])

        f_KAW = f(f > 8e-3 & f <= 2e-2);
        power_spectrum_KAW = power_spectrum(f > 8e-3 & f <= 2e-2);
    end
        
        length_f_MHD = length(f_MHD);
        length_f_KAW = length(f_KAW);
        
        
        if  length_f_MHD > 5 & length_f_KAW > 5
            coeff_MHD = coeffvalues(fit(log(f_MHD.'), log(power_spectrum_MHD.'), 'poly1'));
            coeff_KAW = coeffvalues(fit(log(f_KAW.'), log(power_spectrum_KAW.'), 'poly1')); 
%             coeff_MHD = coeffvalues(fit(f1_slope.', power_spectrum_1_slope.', 'power1')); 
%             coeff_KAW = coeffvalues(fit(f2_slope.', power_spectrum_2_slope.', 'power1'));
        elseif length_f_MHD < 5 & length_f_KAW > 5
            coeff_MHD = [0, 0];
            f_MHD = [];
            power_spectrum_MHD = [];
            coeff_KAW = coeffvalues(fit(log(f_KAW.'), log(power_spectrum_KAW.'), 'poly1'));  
%             coeff_KAW = coeffvalues(fit(f2_slope.', power_spectrum_2_slope.', 'power1'));
        elseif length(f_KAW) < 5 & length_f_MHD > 5
            coeff_MHD = coeffvalues(fit(log(f_MHD.'), log(power_spectrum_MHD.'), 'poly1'));
%             coeff_MHD = coeffvalues(fit(f1_slope.', power_spectrum_1_slope.', 'power1')); 
            coeff_KAW = [0, 0];
            f_KAW = [];
            power_spectrum_KAW = [];
        else
            coeff_MHD = [0, 0];
            f_MHD = [];
            power_spectrum_MHD = [];
            coeff_KAW = [0, 0];
            f_KAW = [];
            power_spectrum_KAW = 0;
        end
 
end