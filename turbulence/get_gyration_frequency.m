function [gyration_frequency] = get_gyration_frequency(magnetometer_data_to_analyze,ion_kg)

    charge = 1.602E-19;
    avagdro_number = 6.022e23;

    gyration_frequency = charge*mean(magnetometer_data_to_analyze)/(ion_kg/avagdro_number);
end