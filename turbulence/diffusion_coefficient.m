function [d] = diffusion_coefficient(power_spectrum_par,f_par)
            d = mean(power_spectrum_par.*f_par); %\delta B_||^2
end
