length_of_window = 600;
max_windows = 10;
slide_length = length_of_window/5;
avg = 10;
mag_is_1_sheath_is_0 = 1;
avoid_boundary_offset = 1;

[location_data] = get_location_data(); 
regions_boundary_data = get_location_regions_boundary_data();
ordered_crossings = crossings_of_interest(regions_boundary_data,mag_is_1_sheath_is_0);
crossing_date = ordered_crossings(8,1);
boundaries_inside = 1;
winds = 0;

fucktuation_perp = zeros(1,50000);
fucktuation_par = zeros(1,50000);
nl_measure = zeros(1,50000);
mva_eigvals = zeros(1,50000);
window_dates = zeros(1,50000);
q = zeros(1,50000);
d = zeros(1,50000);
R = zeros(1,50000);
lat = zeros(1,50000);
LT = zeros(1,50000);
tfb = zeros(1,50000);
current_sheet = zeros(1,50000);

mu_0 = 4*pi()*1e-7; %N/A^2
kB = 1.3806488e-23; %J/K
avagdro_number = 6.022140857e23;    
O_mass = 16/avagdro_number;

for i = 2:169
    start = boundaries_inside;
    [magnetometer_data] = get_magnetometer_data(i);
    [n, length_of_magnetometer_data] = size(magnetometer_data);
    dates = 24*60*(datenum(magnetometer_data(1,:)...
        , magnetometer_data(2,:), magnetometer_data(3,:)...
        , magnetometer_data(4,:), magnetometer_data(5,:)...
        , floor(magnetometer_data(6,:))) - datenum(2004,1,1));
    [boundaries_in_file,boundaries_inside,crossing_date] = find_crossings_in_file(...
        magnetometer_data,dates,ordered_crossings,start,boundaries_inside,crossing_date,mag_is_1_sheath_is_0);       
    [g,number_of_crossings] = size(boundaries_in_file);  

    for j = 1:number_of_crossings 
        [windows] = window_finder(boundaries_in_file,j,number_of_crossings,length_of_window/60,max_windows);
        index_of_crossing = find(dates == boundaries_in_file(8,j));
        if (windows > 1) && ~isempty(index_of_crossing) 
            winds = 5*(min([max_windows,windows])-1)+1;
            sects1 = zeros(winds-avoid_boundary_offset,1);
            sects2 = zeros(winds-avoid_boundary_offset,1);
            sects3 = zeros(winds-avoid_boundary_offset,1);
            sects4 = zeros(winds-avoid_boundary_offset,1);
            sects5 = zeros(winds-avoid_boundary_offset,1);

            sr1 = zeros(winds-avoid_boundary_offset,1);
            slat1 = zeros(winds-avoid_boundary_offset,1);
            sLT1 = zeros(winds-avoid_boundary_offset,1);
            stfb1 = zeros(winds-avoid_boundary_offset,1);
            q1 = zeros(winds-avoid_boundary_offset,1);
            d1 = zeros(winds-avoid_boundary_offset,1);
            cs1 = zeros(winds-avoid_boundary_offset,1);

%            parfor k = avoid_boundary_offset:(winds-1)
           for k = avoid_boundary_offset:(winds-1)
                do_it = false;
                if boundaries_in_file(7,j) == 1 && (index_of_crossing + k*slide_length + length_of_window - 1 <= length_of_magnetometer_data)
                    %go forwards
                    mag_data_to_analyze = magnetometer_data(:, (index_of_crossing + k*slide_length):(index_of_crossing + k*slide_length + length_of_window - 1));
                    do_it =  true;
                elseif boundaries_in_file(7,j) == 2 && (index_of_crossing - k*slide_length - length_of_window + 1 >= 0)
                   %go backwards
                    mag_data_to_analyze = magnetometer_data(:, (index_of_crossing - k*slide_length - length_of_window + 1):(index_of_crossing - k*slide_length));
                    do_it = true; 
                end 

                if do_it
                    loc = location_data(8:11, location_data(1,:) == mag_data_to_analyze(1,1) ...
                        & location_data(2,:) == mag_data_to_analyze(2,1)...
                        & location_data(3,:) == mag_data_to_analyze(3,1) ...
                        & location_data(4,:) == mag_data_to_analyze(4,1));
                    b_r = (1e-9)*mag_data_to_analyze(7,:).';
                    b_theta = (1e-9)*mag_data_to_analyze(8,:).';
                    b_phi = (1e-9)*mag_data_to_analyze(9,:).';

                    %q = no_vitaliy_bullshit_heating_rate_density(b_r,b_theta,b_phi,loc(4),loc(1),1)
                    mag_dates = 24*60*(datenum(mag_data_to_analyze(1,:),mag_data_to_analyze(2,:),...
                        mag_data_to_analyze(3,:),mag_data_to_analyze(4,:),mag_data_to_analyze(5,:),...
                        floor(mag_data_to_analyze(6,:))) - datenum(2004,1,1));

                    B_r = mean(reshape(b_r,[avg,length_of_window/avg])',2);
                    B_theta = mean(reshape(b_theta,[avg,length_of_window/avg])',2);
                    B_phi = mean(reshape(b_phi,[avg,length_of_window/avg,])',2);
                    if ((sum(B_r > 0) > 0 && sum(B_r < 0) > 0) && (sum(B_theta > 0) > 0 && sum(B_theta < 0) > 0)) || ((sum(B_r > 0) > 0 && sum(B_r < 0) > 0) && (sum(B_phi > 0) > 0 && sum(B_phi < 0) > 0)) || ((sum(B_phi > 0) > 0 && sum(B_phi < 0) > 0) && (sum(B_theta > 0) > 0 && sum(B_theta < 0) > 0))
                        cs1(k) = 1;
                    end
                    
%                     [num2str(mag_data_to_analyze(1,1)),' ',num2str(mag_data_to_analyze(2,1)),' ',num2str(mag_data_to_analyze(3,1)),' ',...
%                         num2str(mag_data_to_analyze(4,1)),':',num2str(mag_data_to_analyze(5,1)),':',num2str(mag_data_to_analyze(6,1)),' j =',num2str(j),' k = ',num2str(k)]
% % 
     %                if j == 1428 && k == 12
                    %variance analysis
                    [e1,e2,e3] = B_var(B_r,B_theta,B_phi);
                    v_tot = geomean([e1,e2,e3]);
                    %fluctuation
                    B_vector_mean = get_B_vector_mean(mean(B_r), mean(B_theta), mean(B_phi));
                    [B_fluct_par,B_fluct_perp,B_std_par,B_std_perp] =...
                       get_B_std_vector_components(B_vector_mean, B_r, B_theta, B_phi);
                    fperp_tot = sqrt(mean(sum(B_fluct_perp.^2))); 
                    fpar_tot = sqrt(mean(sum(B_fluct_par.^2)));
                    %nonlinear time series analysis
                    [tot,br1,br2,bt1,bt2,bp1,bp2] = B_nl_embed (B_r,B_theta,B_phi);
                    %heating rate density stuff     
                     H = get_scale_height(loc(4));
                     [ux, uz] = get_v_rel(loc(4));
                     number_density = get_density(loc(1), H, loc(4));                  
                     [q_MHD,dd] = heating_rate_density(B_r,B_theta,B_phi,...
                         ux,uz,number_density,avg,18,H)

                    sects1(k) = mag_dates(1);
                    sects2(k) = tot;
                    sects3(k) = v_tot;
                    sects4(k) = fpar_tot;
                    sects5(k) = fperp_tot;
                    sr1(k) = loc(4); 
                    slat1(k) = loc(1); 
                    sLT1(k) = loc(3); 
                    stfb1(k) = k+1;
                    q1(k) = q_MHD;
                    d1(k) = dd;

%                    if (log10(tot) >= thresh1 || log10(v_tot) >= thresh2 || log10(fpar_tot^2+fperp_tot^2) >= thresh3)
%                       h = figure; set(h,'Visible','off');
%                       set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.7, 0.74]);
%                       plot((0:length(B_theta)-1)*10,B_r,'r','LineWidth',2)
%                       hold on
%                       plot((0:length(B_theta)-1)*10,B_theta,'g','LineWidth',2)
%                       plot((0:length(B_theta)-1)*10,B_phi,'b','LineWidth',2)
%                       ylabel('magnetic field (T)')
%                       xlabel('time (s)')
%                       legend('B_r','B_\theta','B_\phi')
%                       %title([num2str(mag_data_to_analyze(1,1)),'-',num2str(mag_data_to_analyze(2,1)),...
%                       %    '-',num2str(mag_data_to_analyze(3,1)),' ',num2str(mag_data_to_analyze(4,1)),...
%                       %    ':',num2str(mag_data_to_analyze(5,1),'%02.f'),' - NTS: ',...
%                       %    num2str(log10(tot)),' VAR: ',num2str(log10(v_tot)),' FLUC ',num2str(log10(fpar_tot^2+fperp_tot^2))])
%                     title([num2str(mag_data_to_analyze(1,1)),'-',num2str(mag_data_to_analyze(2,1)),...
%                       '-',num2str(mag_data_to_analyze(3,1)),' ',num2str(mag_data_to_analyze(4,1)),...
%                       ':',num2str(mag_data_to_analyze(5,1),'%02.f')])
%                       set(gca, 'FontSize', 24)
%                       set(gca,'LineWidth',2)
%                       saveas(h,[num2str(mag_data_to_analyze(1,1)),'_',num2str(mag_data_to_analyze(2,1)),...
%                           '_',num2str(mag_data_to_analyze(3,1)),'-',num2str(mag_data_to_analyze(4,1)),...
%                           '-',num2str(mag_data_to_analyze(5,1)),'.png']);
% 
%                       h = figure;% set(h,'Visible','off');
%                       set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.7, 0.74]);
%                       subplot(1,3,1)
%                       plot(br1/1e-9,br2/1e-9,'r')
%                       xlabel('B_{r1} (nT)')
%                       ylabel('B_{r2} (nT)')
%                       set(gca, 'FontSize', 24)
%                       set(gca,'LineWidth',2)
%                       subplot(1,3,2)
%                       plot(bt1/1e-9,bt2/1e-9,'g')
%                       xlabel('B_{\theta1} (nT)')
%                       ylabel('B_{\theta2} (nT)')
%                       set(gca, 'FontSize', 24)
%                       set(gca,'LineWidth',2)
%                       subplot(1,3,3)
%                       plot(bp1/1e-9,bp2/1e-9,'b')
%                       xlabel('B_{\phi1} (nT)')
%                       ylabel('B_{\phi2} (nT)')
%                       set(gca, 'FontSize', 24)
%                       set(gca,'LineWidth',2)
%                       saveas(h,['embed_',num2str(mag_data_to_analyze(1,1)),'_',num2str(mag_data_to_analyze(2,1)),...
%                           '_',num2str(mag_data_to_analyze(3,1)),'_',num2str(mag_data_to_analyze(4,1)),...
%                           '_',num2str(mag_data_to_analyze(5,1)),'.png']);
                    %close all
                    %end
                end
            end
            zinds = find(R == 0);

            window_dates(zinds(1):zinds(1)+winds-1-avoid_boundary_offset) = sects1;
            nl_measure(zinds(1):zinds(1)+winds-1-avoid_boundary_offset) = sects2;
            mva_eigvals(zinds(1):zinds(1)+winds-1-avoid_boundary_offset) = sects3;
            fucktuation_par(zinds(1):zinds(1)+winds-1-avoid_boundary_offset) = sects4;
            fucktuation_perp(zinds(1):zinds(1)+winds-1-avoid_boundary_offset) = sects5;
            q(zinds(1):zinds(1)+winds-1-avoid_boundary_offset) = q1;
            d(zinds(1):zinds(1)+winds-1-avoid_boundary_offset) = d1;
            current_sheet(zinds(1):zinds(1)+winds-1-avoid_boundary_offset) = cs1;

            R(zinds(1):zinds(1)+winds-1-avoid_boundary_offset) = sr1;
            lat(zinds(1):zinds(1)+winds-1-avoid_boundary_offset) = slat1;
            LT(zinds(1):zinds(1)+winds-1-avoid_boundary_offset) = sLT1;
            tfb(zinds(1):zinds(1)+winds-1-avoid_boundary_offset) = stfb1;
        end
    end
end
end_cond = R~= 0; 

window_dates = window_dates(end_cond);
nl_measure = nl_measure(end_cond);
fucktuation_perp = fucktuation_perp(end_cond);
fucktuation_par = fucktuation_par(end_cond);
mva_eigvals = mva_eigvals(end_cond);
q = q(end_cond);
d = d(end_cond);
R = R(end_cond);
lat = lat(end_cond);
LT = LT(end_cond);
tfb = tfb(end_cond);
current_sheet = current_sheet(end_cond);

save('window_dates','window_dates')
save('nl_measure','nl_measure')
save('mva_eigvals','mva_eigvals')
save('fucktuation_perp','fucktuation_perp')
save('fucktuation_par','fucktuation_par')
save('q','q')
save('d','d')
save('R','R')
save('lat','lat')
save('LT','LT')
save('tfb','tfb')
save('current_sheet','current_sheet')

%disturbed_anal