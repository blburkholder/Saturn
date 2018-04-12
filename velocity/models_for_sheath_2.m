function [density,temperature,v_r,v_phi] = models_for_sheath_2()
    mag_is_1_sheath_is_2 = 1;

    data = get_LANL_moments();
    boundaries = get_location_regions_boundary_data();
    %sheath from magnetopause
    crossings = crossings_of_interest(boundaries,mag_is_1_sheath_is_2);

    %sheath from shock
    %crossings = crossings_of_interest_2(boundaries);

    resolution_in_minutes = 30;
    slices = (24*60)/resolution_in_minutes;
    k = 60/resolution_in_minutes;
    h = 0;

    save_data = false;
    if save_data
        moment_times_path_Name_w = '/home/computation/Documents/GitProjects/';
        moment_times_file_Name_w = 'sheath_500mins.txt';     
        moment_times_file_Name = horzcat(moment_times_path_Name_w, moment_times_file_Name_w);
        moment_times_fileID = fopen(moment_times_file_Name, 'w');

        header = 'Year        DOY        Hour       Min       Sec        LT\n'; 
        fprintf(moment_times_fileID, header); 
    end

    dates = 24*60*(datenum(data(2,:),1,1) + (data(3,:)-1) + data(4,:)/24 + data(5,:)/(24*60)...
            + data(6,:)/(24*60*60) - datenum(2004,1,1));

    [x,y] = size(crossings);
    crossing_by_crossing_average = zeros(y,5);
    cbcLT = zeros(y,1);
    %cbcLT2 = zeros(y,1);
    cbcR = zeros(y,1);
    %cbc = zeros(y,1);

    figure
    hold on

    for i = 1:y
        %can keep it this way even if you choose crossings_of_interest2
        if mag_is_1_sheath_is_2 == 2
            if crossings(7,i) == 2
                ze_condition = ~isnan(data(8,:)) & dates >= crossings(8,i) &...
                dates <= crossings(8,i) + crossings(9,i)/2 & abs(data(27,:)) < 30 & ~data(29,:) &...
                ~data(30,:) & ((data(28,:) < 12 & data(37,:)) | (data(28,:) >= 12 & ~data(37,:))) & dates - crossings(8,i) <= 100;

                %corotation viewing flag comparison
                %ze_condition2 = ~isnan(data(8,:)) & dates >= crossings(8,i) &...
                %dates <= crossings(8,i) + crossings(9,i)/2 & abs(data(27,:)) < 30 & ~data(29,:) &...
                %~data(30,:) & dates - crossings(8,i) <= 500;
            else
                ze_condition = ~isnan(data(8,:)) & dates <= crossings(8,i) &...
                dates >= crossings(8,i) - crossings(9,i)/2 & abs(data(27,:)) < 30 & ~data(29,:) &...
                ~data(30,:) & ((data(28,:) < 12 & data(37,:)) | (data(28,:) >= 12 & ~data(37,:))) & crossings(8,i) - dates <= 100;

                %ze_condition2 = ~isnan(data(8,:)) & dates <= crossings(8,i) &...
                %dates >= crossings(8,i) - crossings(9,i)/2 & abs(data(27,:)) < 30 & ~data(29,:) &...
                %~data(30,:) & crossings(8,i) - dates <= 500;
            end
        else
            if crossings(7,i) == 1
                ze_condition = ~isnan(data(8,:)) & dates >= crossings(8,i) &...
                dates <= crossings(8,i) + crossings(9,i)/2 & abs(data(27,:)) < 30 & ~data(29,:) &...
                ~data(30,:) & ~data(37,:) & abs(dates - crossings(8,i)) <= 100;
            else
                ze_condition = ~isnan(data(8,:)) & dates <= crossings(8,i) &...
                dates >= crossings(8,i) - crossings(9,i)/2 & abs(data(27,:)) < 30 & ~data(29,:) &...
                ~data(30,:) & ~data(37,:) & abs(crossings(8,i) - dates) <= 100;
            end
        end

        densities = data(7,ze_condition & data(7,:) > 0);
        temps = data(8,ze_condition & data(8,:) > -999);
        v_r_rel = data(9,ze_condition & data(9,:) ~= -999);
        v_r_fluct = sqrt((v_r_rel - mean(v_r_rel)).^2);
        v_phi_rel = data(10,ze_condition & data(10,:) ~= -999);
        %v_phi_rel2 = data(10,ze_condition2 & data(10,:) ~= -999);
        
        LT = data(28, ze_condition & data(10,:) ~= -999);
        %LT2 = data(28, ze_condition & data(10,:) ~= -999);
        R = data(26, ze_condition & data(10,:) ~= -999);

        %avoid repeat data
        data(7,ze_condition) = 0;
        data(8,ze_condition) = -999;
        data(9,ze_condition) = -999;
        data(10,ze_condition) = -999;

        %calculates an average for each boundary
        if ~isempty(v_phi_rel)
            %cbc(i) = nanmean(v_phi_rel2);
            %cbcLT2(i) = nanmean(LT2);
            crossing_by_crossing_average(i,1) = nanmean(v_phi_rel);
            crossing_by_crossing_average(i,2) = nanmean(v_r_rel);
            crossing_by_crossing_average(i,3) = nanmean(v_r_fluct);
            crossing_by_crossing_average(i,4) = geomean(densities);
            crossing_by_crossing_average(i,5) = nanmean(temps);
            if crossings(7,i) == 2
                cbcLT(i) = LT(1);
                cbcR(i) = R(1);
            else
                cbcLT(i) = LT(end);
                cbcR(i) = R(end);
            end
            
            if mag_is_1_sheath_is_2
                [normalizing_vphi,~] = get_v_rel(R);
            else
                normalizing_vphi = ones(size(v_phi_rel));
            end

            dynpress = get_dynamic_pressure(cbcR(i),pi*cbcLT(i)/12 - pi);
            if dynpress < 0.01
                scatter(ones(size(v_phi_rel))*cbcLT(i),1000*v_phi_rel./normalizing_vphi,'r.');
            elseif dynpress >= 0.01 && dynpress < 0.05
                scatter(ones(size(v_phi_rel))*cbcLT(i),1000*v_phi_rel./normalizing_vphi,'g.');
            else
                scatter(ones(size(v_phi_rel))*cbcLT(i),1000*v_phi_rel./normalizing_vphi,'b.');  
            end    
        end
  
        if save_data
            these_dates = data(2:6,ze_condition & data(9,:) ~= -999);
            pre_and_post_noon = these_dates(:,LT > 10 & LT < 14);
            pre_and_post_noon_times = LT(LT > 10 & LT < 14);
            h = h  + length(pre_and_post_noon_times);
            for not_k = 1:length(pre_and_post_noon_times)
                times_to_write = [pre_and_post_noon(:,not_k)',pre_and_post_noon_times(not_k)]; 
                formatSpec_w = ['%f %f %f %f %f %f\n'];          
                fprintf(moment_times_fileID, formatSpec_w, times_to_write);
            end
        end
    end

%     cbc_dynpress = zeros(y,1);
%     for g = 1:y
%         cbc_dynpress(g) = get_dynamic_pressure(cbcR(g),pi*cbcLT(g)/12 - pi);
%     end

%     figure
%     scatter(cbcLT(cbc_dynpress < 0.01),crossing_by_crossing_average(cbc_dynpress < 0.01,1),'r.')
%     hold on
%     scatter(cbcLT(cbc_dynpress >= 0.01 & cbc_dynpress < 0.05),crossing_by_crossing_average(cbc_dynpress >= 0.01 & cbc_dynpress < 0.05 ,1),'g.')
%     scatter(cbcLT(cbc_dynpress >= 0.05),crossing_by_crossing_average(cbc_dynpress >= 0.05,1),'b.')


    %flag_view = zeros(1,length(all_v_r(all_vr_LT < 14 & all_vr_LT > 10)));
    %flag_edge_anode = zeros(1,length(all_v_r(all_vr_LT < 14 & all_vr_LT > 10)));
    %fid = fopen('/home/computation/GitProjects/pre&postnoon_fixed&flagged&halved.txt');
    %FC = textscan(fid,'%s','Delimiter','\n');
    %for i = 1:length(all_v_r(all_vr_LT < 14 & all_vr_LT > 10))
    %    b = FC{1}(i);
    %    flag_view(i) = b{1}(1);
    %    flag_edge_anode(i) = b{1}(2);
    %end

    %all_v_r = all_v_r((flag_view == 105 | flag_view == 101));
    %all_vr_fluct = all_vr_fluct((flag_view == 105 | flag_view == 101));
    %all_vr_LT = all_vr_LT((flag_view == 105 | flag_view == 101));
    %all_v_phi = all_v_phi((flag_view == 105 | flag_view == 101));
    %all_vphi_LT = all_vphi_LT((flag_view == 105 | flag_view == 101));
      
    avg_model_data_T = nan(slices,1);
    avg_model_data_density = nan(slices,1);
    avg_model_data_v_r = nan(slices,1);
    avg_model_data_v_phi = nan(slices,1);
    %avg_model_data_v_phi2 = nan(slices,1);
    vr_fluct = zeros(slices,1);

    for i = 1:slices
        avg_model_data_density(i) = geomean(crossing_by_crossing_average(floor(k*cbcLT) == i & crossing_by_crossing_average(:,4),4));
        avg_model_data_T(i) =  mean(crossing_by_crossing_average(floor(k*cbcLT) == i & crossing_by_crossing_average(:,5),5));
        avg_model_data_v_r(i) =  mean(crossing_by_crossing_average(floor(k*cbcLT) == i & crossing_by_crossing_average(:,2),2));
        avg_model_data_v_phi(i) = mean(crossing_by_crossing_average(floor(k*cbcLT) == i & crossing_by_crossing_average(:,1),1));
        %avg_model_data_v_phi2(i) = mean(cbc(floor(k*cbcLT2) == i & cbc));
        vr_fluct(i) = mean(crossing_by_crossing_average(floor(k*cbcLT) == i,3));
    end

%    avg_model_data_T = nan(50,1);
%    avg_model_data_density = nan(50,1);
%    avg_model_data_v_r = nan(50,1);
%    avg_model_data_v_phi = nan(50,1);
%    vr_fluct = zeros(50,1);

%    for i = 1:100
%        avg_model_data_density(i) = geomean(crossing_by_crossing_average(floor(k*cbcR) == i & crossing_by_crossing_average(:,4),4));
%        avg_model_data_T(i) =  mean(crossing_by_crossing_average(floor(k*cbcR) == i & crossing_by_crossing_average(:,5),5));
%        avg_model_data_v_r(i) =  mean(crossing_by_crossing_average(floor(k*cbcR) == i & crossing_by_crossing_average(:,2),2));
%        avg_model_data_v_phi(i) = mean(crossing_by_crossing_average(floor(k*cbcR) == i & crossing_by_crossing_average(:,1),1));
%        vr_fluct(i) = mean(crossing_by_crossing_average(floor(k*cbcR) == i,3));
%    end

    density_first = find(~isnan(avg_model_data_density),1,'first');
    density_last = find(~isnan(avg_model_data_density),1,'last');
    temp_first = find(~isnan(avg_model_data_T),1,'first');
    temp_last = find(~isnan(avg_model_data_T),1,'last');
    v_r_first = find(~isnan(avg_model_data_v_r),1,'first');
    v_r_last = find(~isnan(avg_model_data_v_r),1,'last');
    v_phi_first = find(~isnan(avg_model_data_v_phi),1,'first');
    v_phi_last = find(~isnan(avg_model_data_v_phi),1,'last');
    %v_phi_first2 = find(~isnan(avg_model_data_v_phi2),1,'first');
    %v_phi_last2 = find(~isnan(avg_model_data_v_phi2),1,'last');

    interp_density = interp1(find(~isnan(avg_model_data_density)),avg_model_data_density(~isnan(avg_model_data_density)),density_first:density_last);
    interp_temp = interp1(find(~isnan(avg_model_data_T)),avg_model_data_T(~isnan(avg_model_data_T)),temp_first:temp_last);
    interp_vr = interp1(find(~isnan(avg_model_data_v_r)),avg_model_data_v_r(~isnan(avg_model_data_v_r)),v_r_first:v_r_last);
    interp_vphi = interp1(find(~isnan(avg_model_data_v_phi)),avg_model_data_v_phi(~isnan(avg_model_data_v_phi)),v_phi_first:v_phi_last);
    
    smooth_density = transpose(smooth(interp_density));
    smooth_temp = transpose(smooth(interp_temp));
    smooth_vr = transpose(smooth(interp_vr));
    smooth_vphi = transpose(smooth(interp_vphi));

    %smooth_density = interp_density;
    %smooth_temp = interp_temp;
    %smooth_vr = interp_vr;
    %smooth_vphi = interp_vphi;

    p1 = polyfit(temp_first:temp_last,smooth_temp,2);
    p2 = polyfit(density_first:density_last,smooth_density,2);
    p22 = polyfit(density_first+4:density_last-6,smooth_density(4:end-7),2);
    p3 = polyfit(v_r_first:v_r_last,smooth_vr,2);
    p4 = polyfit(v_phi_first:v_phi_last,smooth_vphi,3);
    p44 = polyfit([7 8 9 10 24 38 39 40 41],[-200 -200 -200 -200 0 200 200 200 200],3);

avg_model_data_v_phi(v_phi_first:v_phi_last) = interp_vphi;

    d2 = avg_model_data_v_phi(v_phi_first:22);
    %d22 = avg_model_data_v_phi2(v_phi_first2+4:22);

%LOTS OF PLOTS
%---------------------------------------------------------------------------
 %   figure
 %   h1 =scatter(1:14,d2,'b')
 %   set(h1, 'SizeData', 100)
 %   hold on
 %   h2 = scatter(1:14,d22,'r*')
 %   set(h2, 'SizeData', 100)
 %   xlabel('Local Time','Interpreter','latex')
 %   ylabel('km/s','Interpreter','latex')
 %   legend('anticorotation assumption','viewing ignored')

    d3 = avg_model_data_v_phi(26:v_phi_last);
    figure
    bar([slices/4+2,3*slices/4-2],[mean(d2),mean(d3)],0.8)
    hold on
    %h = boxplot([d2(al+1:end),d3],'positions',[15,33],'Whisker',0,'Colors','k');
    %set(h,{'linew'},{1.5});
    %set(h(7,:),'Visible','off');
    xlim auto
    title('v_{\phi}');
    line([slices/2 slices/2],[-300,300],'Color','k');
    line([0 slices],[0,0],'Color','k');
    %plot(v_phi_first:v_phi_last,smooth_vphi)
    %scatter(v_phi_first:v_phi_last,avg_model_data_v_phi(v_phi_first:v_phi_last),[],vr_fluct(v_phi_first:v_phi_last),'.','SizeData',500)
    scatter(v_phi_first:v_phi_last,avg_model_data_v_phi(v_phi_first:v_phi_last),'.','SizeData',500)
    %a = plot(v_phi_first:23,-avg_model_data_v_phi(v_phi_first:23),'r--');
    %b = plot(48 - v_phi_last:23,flipud(avg_model_data_v_phi(25:v_phi_last)),'b--');
    c = plot(v_phi_first:23,avg_model_data_v_phi(v_phi_first:23),'r--');
    d = plot(25:v_phi_last,avg_model_data_v_phi(25:v_phi_last),'b--');
    ylabel('km/s')
    xlabel('Local Time')
    %h = colorbar;
    %ylabel(h,'v_r fluctuation')
    ax = gca;
    ax.XTick = [0 6 12 18 24 30 36 42];
    ax.XTickLabel = {'0',' ','6:00',' ','12:00',' ','18:00',' '};
    %x1 = linspace(v_phi_first,v_phi_last);
    %y2 = polyval(p4,x1);
    %plot(x1,y2)
    %y2 = polyval(p44,x1);
    %plot(x1,y2,'g')
   
    %rgb_colors = parula(100);
    %r = rgb_colors(:,1);
    %g = rgb_colors(:,2);
    %b = rgb_colors(:,3);

    %bb2 = bar(15,nanmean(d2),14);
    %alpha(bb2,0.1);
    %set(bb2(1),'facecolor',[r(ceil(nanmean(vr_fluct(v_phi_first:22))/maximum*100)) g(ceil(nanmean(vr_fluct(v_phi_first:22))/maximum*100)) b(ceil(nanmean(vr_fluct(v_phi_first:22))/maximum*100))]); 
    %bb3 = bar(33,nanmean(d3),14);
    %alpha(bb3,0.1);
    %set(bb3(1),'facecolor',[r(ceil(nanmean(vr_fluct(26:v_phi_last))/maximum*100)) g(ceil(nanmean(vr_fluct(26:v_phi_last))/maximum*100)) b(ceil(nanmean(vr_fluct(26:v_phi_last))/maximum*100))]); 
%{
    figure
    %plot(v_r_first:v_r_last,smooth_vr)
    scatter(v_r_first:v_r_last,avg_model_data_v_r(v_r_first:v_r_last),'.','SizeData',500)
    line([slices/2 slices/2],[-50,250],'Color','k')
    line([0 slices],[0,0],'Color','k');
    title('v_r')
    hold on
    x1 = linspace(v_r_first,v_r_last);
    y1 = polyval(p3,x1);
    plot(x1,y1)
    x1 = v_r_first:v_r_last;
    a = 0.5*p3(2)/p3(1)+24;
    y2 = p3(1)*x1.^2 + (p3(2) - 2*p3(1)*a)*x1 + p3(1)*a^2-p3(2)*a+p3(3);
    p3(2) = p3(2) - 2*p3(1)*a;
    p3(3) = p3(1)*a^2-p3(2)*a+p3(3);
    plot(x1,y2,'g')
    ax = gca;
    ax.XTick = [0 10 20 30 40 50];
    ax.XTickLabel = {'0','5:00','10:00','15:00','20:00','-'};
    hold off

    figure
    %plot(density_first:density_last,smooth_density)
    scatter(density_first:density_last,avg_model_data_density(density_first:density_last),'.','SizeData',500)
    line([slices/2 slices/2],[0,0.3],'Color','k')
    title('density')
    hold on
    x1 = linspace(density_first,density_last);
    y2 = polyval(p2,x1);
    plot(x1,y2)
    x1 = density_first:density_last;
    a = 0.5*p22(2)/p22(1)+24;
    y2 = p22(1)*x1.^2 + (p22(2) - 2*p22(1)*a)*x1 + p22(1)*a^2-p22(2)*a+p22(3);
    p22(2) = p22(2) - 2*p22(1)*a;
    p22(3) = p22(1)*a^2-p22(2)*a+p22(3);
    plot(x1,y2,'g')
    ax = gca;
    ax.XTick = [0 10 20 30 40 50];
    ax.XTickLabel = {'0','5:00','10:00','15:00','20:00','-'};
    hold off

    figure
    %plot(1:72,avg_model_data_T(1:72))
    %plot(temp_first:temp_last,smooth_temp)
    scatter(temp_first:temp_last,avg_model_data_T(temp_first:temp_last),'.','SizeData',500)
    title('temperature')
    line([slices/2 slices/2],[0,600],'Color','k')
    hold on
    x1 = linspace(temp_first,temp_last);
    y2 = polyval(p1,x1);
    plot(x1,y2)
    x1 = temp_first:temp_last;
    a = 0.5*p1(2)/p1(1)+24;
    y2 = p1(1)*x1.^2 + (p1(2) - 2*p1(1)*a)*x1 + p1(1)*a^2-p1(2)*a+p1(3);
    p1(2) = p1(2) - 2*p1(1)*a;
    p1(3) = p1(1)*a^2-p1(2)*a+p1(3);
    plot(x1,y2,'g')
    ax = gca;
    ax.XTick = [0 10 20 30 40 50];
    ax.XTickLabel = {'0','5:00','10:00','15:00','20:00','-'};
    hold off
%}
    %figure
    %plot(interp_density.*interp_vphi);
    %line([slices/2 slices/2],[-40,100],'Color','r')
    %title('Mass flux v_{phi}*rho')
    %line([0 40],[0,0],'Color','r')
%--------------------------------------------------------------------------

    temperature = p1;
    density = p22;
    v_r = p3;
    v_phi = p44;