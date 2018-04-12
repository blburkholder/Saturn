data = get_LANL_moments();
boundaries = get_location_regions_boundary_data();
mag_is_1_sheath_is_2 = 1;
crossings = crossings_of_interest(boundaries,mag_is_1_sheath_is_2);
dates = 24*60*(datenum(data(2,:),1,1) + (data(3,:)-1) + data(4,:)/24 + data(5,:)/(24*60)...
            + data(6,:)/(24*60*60) - datenum(2004,1,1));
[x,y] = size(crossings);
avg = zeros(113,35,4);
dawn_0_dusk_1 = zeros(71,1);
%avg = zeros(79,35,4);
%dawn_0_dusk_1 = zeros(79,1);
j = 0;
    for i = 1:y
        if mag_is_1_sheath_is_2 == 1
            if crossings(7,i) == mag_is_1_sheath_is_2
                ze_condition = ~isnan(data(8,:)) & dates >= crossings(8,i) &...
                dates <= crossings(8,i) + crossings(9,i)/2 & abs(data(27,:)) < 30 & ~data(29,:) &...
                ~data(30,:) & ~data(37,:) & dates - crossings(8,i) <= 200;
            else
                ze_condition = ~isnan(data(8,:)) & dates <= crossings(8,i) &...
                dates >= crossings(8,i) - crossings(9,i)/2 & abs(data(27,:)) < 30 & ~data(29,:) &...
                ~data(30,:) & ~data(37,:) & crossings(8,i) - dates <= 200;
            end
        else
            if crossings(7,i) == mag_is_1_sheath_is_2 == 2
%not implemented            
            else
            
            end
        end

        densities = data(7,ze_condition & data(7,:) > 0);
        temps = data(8,ze_condition & data(8,:) > -999);
        %v_r_rel = data(9,ze_condition & data(9,:) ~= -999);
        %v_phi_rel = data(10,ze_condition & data(10,:) ~= -999);
        LT = data(28, ze_condition & data(10,:) ~= -999);
        times  = dates(ze_condition & data(7,:) > 0);

        if length(densities) > 20
            times = times - times(1);    
            xq = 0:2:2*floor(times(end)/2);
            df = 1/(length(xq)*2);
            f = [0:length(xq)/2-1]*df;

            %subtract off DC component
            entropy = densities.^(-2/3).*temps;
            entropy = entropy - mean(entropy);
            densities = densities - mean(densities);
            temps =  temps - mean(temps);

            vq1 = interp1(times,entropy,xq);
            vq2 = interp1(times,temps,xq);
            vq3 = interp1(times,densities,xq);
            %vq = interp1(times,v_r_rel,xq);
            %vq = interp1(times,v_phi_rel,xq);

            X = fft(vq1);
            Y = fft(vq2);
            Z = fft(vq3);
            j = j+1;
            avg(j,:,1) = abs(X(1:35));
            avg(j,:,2) = abs(Y(1:35));
            avg(j,:,3) = abs(Z(1:35));
            avg(j,:,4) = f(1:35);
            if mean(LT >=12)
                dawn_0_dusk_1(j) = 1;
            end
            %figure
            %semilogy(abs(Y))

            %figure
            %plot(xq,vq);
            %hold on
            %plot(times,v_r_rel);
        end
    end

freq_avg = nan(j,20,3);
for k = 1:j
    for i = 1:20
        freq_avg(k,i,1) = nanmean(avg(k,avg(k,:,4) >= (i-1)*0.01 & avg(k,:,4) < i*0.01,1));
        freq_avg(k,i,2) = nanmean(avg(k,avg(k,:,4) >= (i-1)*0.01 & avg(k,:,4) < i*0.01,2));
        freq_avg(k,i,3) = nanmean(avg(k,avg(k,:,4) >= (i-1)*0.01 & avg(k,:,4) < i*0.01,3));
    end
end


figure
plot(nanmean(freq_avg(dawn_0_dusk_1 == 1,:,1)))
hold on
plot(nanmean(freq_avg(dawn_0_dusk_1 == 0,:,1)))
title('entropy')
legend('dusk','dawn');

figure
plot(nanmean(freq_avg(dawn_0_dusk_1 == 1,:,2)))
hold on
plot(nanmean(freq_avg(dawn_0_dusk_1 == 0,:,2)))
title('temperature')
legend('dusk','dawn')

figure
plot(nanmean(freq_avg(dawn_0_dusk_1 == 1,:,3)))
hold on
plot(nanmean(freq_avg(dawn_0_dusk_1 == 0,:,3)))
title('density')
legend('dusk','dawn')


