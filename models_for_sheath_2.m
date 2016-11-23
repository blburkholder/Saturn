function [density,temperature,v_r,v_phi] = models_for_sheath_2()  
    data = get_LANL_moments();
    boundaries = get_location_regions_boundary_data();
    crossings = crossings_of_interest(boundaries,0);
    %crossings = crossings_of_interest_2(boundaries);
    %crossings = horzcat(crossings1,crossings2);
    resolution_in_minutes = 30;
    slices = (24*60)/resolution_in_minutes;
    k = 60/resolution_in_minutes; 

    dates = 24*60*(datenum(data(2,:),1,1) + (data(3,:)-1) + data(4,:)/24 + data(5,:)/(24*60)...
            + data(6,:)/(24*60*60) - datenum(2004,1,1));

    [x,y] = size(crossings);
    crossing_by_crossing_average = zeros(y,5);
    cbcLT = zeros(y,1);
    for i = 1:y
        if crossings(7,i) == 2
            ze_condition = ~isnan(data(8,:)) & dates >= crossings(8,i) & dates <= crossings(8,i) + crossings(9,i) & abs(data(27,:)) < 30 & ~data(29,:) & ~data(30,:) & ((data(28,:) < 12 & data(37,:)) | (data(28,:) >= 12 & ~data(37,:))) & dates - crossings(8,i) <= 500;
        else
            ze_condition = ~isnan(data(8,:)) & dates <= crossings(8,i) & dates >= crossings(8,i) - crossings(9,i) & abs(data(27,:)) < 30 & ~data(29,:) & ~data(30,:) & ((data(28,:) < 12 & data(37,:)) | (data(28,:) >= 12 & ~data(37,:))) & crossings(8,i) - dates <= 500;
        end

        densities = data(7,ze_condition & data(7,:) > 0);
        temps = data(8,ze_condition & data(8,:) > -999);
        v_r_rel = data(9,ze_condition & ~isnan(data(8,:)) & data(9,:) ~= -999);
        v_r_fluct = sqrt((v_r_rel - mean(v_r_rel)).^2);
        v_phi_rel = data(10,ze_condition & ~isnan(data(8,:)) & data(10,:) ~= -999);
        LT = data(28, ze_condition & ~isnan(data(8,:)) & data(10,:) ~= -999);

        %avoid repeat data
        data(7,ze_condition) = 0;
        data(8,ze_condition) = -999;
        data(9,ze_condition) = -999;
        data(10,ze_condition) = -999;

        if ~isempty(v_phi_rel)
            crossing_by_crossing_average(i,1) = nanmean(v_phi_rel);
            crossing_by_crossing_average(i,2) = nanmean(v_r_rel);
            crossing_by_crossing_average(i,3) = nanmean(v_r_fluct);
            crossing_by_crossing_average(i,4) = geomean(densities);
            crossing_by_crossing_average(i,5) = nanmean(temps);
            cbcLT(i) = nanmean(LT);
        end  
    end

    %flag_view = zeros(1,length(all_v_r(all_vr_LT < 14 & all_vr_LT > 10)));
    %flag_edge_anode = zeros(1,length(all_v_r(all_vr_LT < 14 & all_vr_LT > 10)));
    %fid = fopen('/home/computation/GitProjects/pre&postnoon_fixed&flagged.txt');
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
    vr_fluct = zeros(slices,1);

    for i = 1:slices
        avg_model_data_density(i) = geomean(crossing_by_crossing_average(floor(k*cbcLT) == i & crossing_by_crossing_average(:,4),4));
        avg_model_data_T(i) =  mean(crossing_by_crossing_average(floor(k*cbcLT) == i & crossing_by_crossing_average(:,5),5));
        avg_model_data_v_r(i) =  mean(crossing_by_crossing_average(floor(k*cbcLT) == i & crossing_by_crossing_average(:,2),2));
        avg_model_data_v_phi(i) = mean(crossing_by_crossing_average(floor(k*cbcLT) == i & crossing_by_crossing_average(:,1),1));
        vr_fluct(i) = mean(crossing_by_crossing_average(floor(k*cbcLT) == i,3));
    end

    density_first = find(~isnan(avg_model_data_density),1,'first');
    density_last = find(~isnan(avg_model_data_density),1,'last');
    temp_first = find(~isnan(avg_model_data_T),1,'first');
    temp_last = find(~isnan(avg_model_data_T),1,'last');
    v_r_first = find(~isnan(avg_model_data_v_r),1,'first');
    v_r_last = find(~isnan(avg_model_data_v_r),1,'last');
    v_phi_first = find(~isnan(avg_model_data_v_phi),1,'first');
    v_phi_last = find(~isnan(avg_model_data_v_phi),1,'last');

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


    d2 = avg_model_data_v_phi(v_phi_first+4:22);
    d3 = avg_model_data_v_phi(26:v_phi_last);

    figure
    %h = boxplot([d2,d3],'positions',[15,33],'Whisker',0,'Colors','k');
    %set(h,{'linew'},{1.5});
    %set(h(7,:),'Visible','off');
    hold on
    xlim auto
    title('v_{\phi}');
    line([slices/2 slices/2],[-300,300],'Color','k');
    line([0 slices],[0,0],'Color','k');
    %scatter(v_phi_first:v_phi_last,avg_model_data_v_phi(v_phi_first:v_phi_last),[],vr_fluct(v_phi_first:v_phi_last),'.','SizeData',500)
    scatter(v_phi_first+4:v_phi_last,avg_model_data_v_phi(v_phi_first+4:v_phi_last),'.','SizeData',500)
    ylabel('km/s')
    xlabel('Local Time')
    %h = colorbar;
    %ylabel(h,'v_r fluctuation')
    ax = gca;
    ax.XTick = [0 10 20 30 40 50];
    ax.XTickLabel = {'0','5:00','10:00','15:00','20:00','-'};
    saveas(gcf,'\home\computation\Pictures\total_vphi_scatter_flags','jpg');
    x1 = linspace(v_phi_first,v_phi_last);
    y2 = polyval(p4,x1);
    plot(x1,y2)
    y2 = polyval(p44,x1);
    plot(x1,y2,'g')

   
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

    figure
    plot(v_r_first:v_r_last,smooth_vr)
    line([slices/2 slices/2],[0,200],'Color','k')
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
    plot(density_first:density_last,smooth_density)
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
    plot(temp_first:temp_last,smooth_temp)
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

    %figure
    %plot(interp_density.*interp_vphi);
    %line([slices/2 slices/2],[-40,100],'Color','r')
    %title('Mass flux v_{phi}*rho')
    %line([0 40],[0,0],'Color','r')

    temperature = p1;
    density = p22;
    v_r = p3;
    v_phi = p44;