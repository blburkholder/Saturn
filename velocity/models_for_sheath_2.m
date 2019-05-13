    mag_is_1_sheath_is_2 = 1;
    avagdro_number = 6.022140857e23;
    mu_0 = 4*pi()*1e-7; %N/A^2
    kB = 1.3806488e-23; %J/K


    load(['extendo_moms.mat'])
    load(['q.mat'])
    load(['d.mat'])
    load(['window_dates.mat'])
    data = extendo_moms;
    boundaries = get_location_regions_boundary_data();
    %sheath from magnetopause
    crossings = crossings_of_interest(boundaries,mag_is_1_sheath_is_2);

    %sheath from shock
    %crossings = crossings_of_interest_2(boundaries);

    resolution_in_minutes = 30;
    slices = (24*60)/resolution_in_minutes;
    k = 60/resolution_in_minutes;

    dates = 24*60*(datenum(data(2,:),1,1) + (data(3,:)-1) + data(4,:)/24 + data(5,:)/(24*60)...
            + data(6,:)/(24*60*60) - datenum(2004,1,1));

    [x,y] = size(crossings);
    crossing_by_crossing_average = zeros(y,4);
    cbcLT = zeros(y,1);
    %cbcR = zeros(y,1);

    grr_LT = [];
    grr_q = [];
    grr_d = [];
    grr_vphi = [];

    for i = 1:y
        %can keep it this way even if you choose crossings_of_interest2
        if mag_is_1_sheath_is_2 == 2
            if crossings(7,i) == 2
                ze_condition = dates >= crossings(8,i) & data(end-2,:) ~= -5 &...
                    dates <= crossings(8,i) + crossings(9,i)/2 & abs(data(27,:)) < 30 &...
                    ((data(28,:) < 12 & data(37,:)) | (data(28,:) >= 12 & ~data(37,:))) &...
                    dates - crossings(8,i) <= 100;
            else
                ze_condition = dates <= crossings(8,i) & data(end-2,:) ~= -5 &...
                    dates >= crossings(8,i) - crossings(9,i)/2 & abs(data(27,:)) < 30 &...
                    ((data(28,:) < 12 & data(37,:)) | (data(28,:) >= 12 & ~data(37,:))) &...
                    crossings(8,i) - dates <= 100;
            end
        else
            if crossings(7,i) == 1
                ze_condition = dates >= crossings(8,i) & data(end-2,:) ~= -5 &...
                    dates <= crossings(8,i) + crossings(9,i)/2 & abs(data(27,:)) < 30 &...
                    ~data(37,:) & abs(dates - crossings(8,i)) <= 100;
            else
                ze_condition = dates <= crossings(8,i) & data(end-2,:) ~= -5 &...
                    dates >= crossings(8,i) - crossings(9,i)/2 & abs(data(27,:)) < 30 &...
                    ~data(37,:) & abs(crossings(8,i) - dates) <= 100;
            end
        end
        if mag_is_1_sheath_is_2 == 1
            ze_condition2 = data(10,:) ~= -999 & data(12,:) > 0 & data(13,:) ~= -999;
            ze_condition22 = data(20,:) ~= -999 & data(22,:) > 0  & data(23,:) ~= -999;
            %num_density = (1e6)*(data(12,ze_condition & ze_condition2));
            %mass_density = num_density*O_mass;
            v_phi_rel = data(10,ze_condition & ze_condition2);
            v_phi_relO = data(20,ze_condition & ze_condition22);
            temps = data(13,ze_condition & ze_condition2);
        else
            ze_condition2 = data(10,:) ~= -999 & data(7,:) > 0 & data(8,:) ~= -999;
            ze_condition22 = data(20,:) ~= -999 & data(17,:) > 0& data(18,:) ~= -999;
            %num_density = (1e6)*(data(7,ze_condition &ze_condition2));
            %mass_density = num_density*H_mass;
            v_phi_rel = data(10,ze_condition & ze_condition2);
            v_phi_relO = data(20,ze_condition & ze_condition22);
            temps = data(8,ze_condition & ze_condition2);
        end
        bx = data(end-2,ze_condition & ze_condition2);
        by = data(end-3,ze_condition & ze_condition2);
        bz = data(end-4,ze_condition & ze_condition2);
        LT = data(28, ze_condition & ze_condition2);
        %R = data(26, ze_condition & data(10,:) ~= -999);

        %avoid repeat data
        data(7,ze_condition) = 0;
        data(8,ze_condition) = -999;
        data(10,ze_condition) = -999;
        data(12,ze_condition) = 0;
        data(13,ze_condition) = -999;
        data(20,ze_condition) = -999;

       % if crossings(7,i) == 1
       %     z_condition = window_dates >= crossings(8,i) & window_dates <= crossings(8,i) + crossings(9,i)/2;
       % else
       %     z_condition = window_dates <= crossings(8,i) & window_dates >= crossings(8,i) - crossings(9,i)/2;
       % end

        %calculates an average for each boundary
        if ~isempty(v_phi_rel)% && sum(z_condition) > 0
            bsqr = bx(:).^2+by(:).^2+bz(:).^2;

            %qq = q(z_condition);
            %dd = d(z_condition);
            %window_dates(z_condition)

            %crossing_by_crossing_average(i,1) = max(q(z_condition));
            %crossing_by_crossing_average(i,2) = max(d(z_condition));
            crossing_by_crossing_average(i,3) = nanmean(v_phi_rel);
            crossing_by_crossing_average(i,4) = nanmean(v_phi_relO);
            if crossings(7,i) == 2
                cbcLT(i) = LT(1);
             %   cbcR(i) = R(1);
            else
                cbcLT(i) = LT(end);
             %   cbcR(i) = R(end);
            end
            %scatter(LT,v_phi_rel_O,[],log10(q),'.');
            %scatter(cbcLT(i),mean(v_phi_rel),'b.');
%             if cbcLT(i) > 12
%                 crossing_by_crossing_average(i,3) = nanmean(v_phi_rel) - 200;
%             else
%                 crossing_by_crossing_average(i,3) = nanmean(v_phi_rel) + 200;
%             end
        end
i
    end

% figure
% q = crossing_by_crossing_average(:,1);
% d = crossing_by_crossing_average(:,2);
% scatter(d,q,'k+')
% xlabel('log_{10}D_\perp (m^2/s)')
% ylabel('log_{10}q (W/m^3)')
% set(gca, 'FontSize', 16)


figure
vphi = crossing_by_crossing_average(:,3);
vphiO = crossing_by_crossing_average(:,4);
plot([3 20],[0 0],'k')
hold on
plot([12 12],[-400 500],'k')
p1 = scatter(cbcLT(vphi<1000),vphi(vphi<1000),'g.');
p2 = scatter(cbcLT(vphiO<1000),vphiO(vphiO<1000),'b.');
legend([p1,p2],'protons','W^+')
if mag_is_1_sheath_is_2 == 1
    scatter(cbcLT(vphi<=0),vphi(vphi<=0),'rp')
    scatter(cbcLT(vphiO<=0),vphiO(vphiO<=0),'rp')
    title('Averaged Magnetospheric Flows')
    axis([3 20 -200 500])
else
    scatter(cbcLT((vphi<=0 & cbcLT>12)|(vphi>=0 & cbcLT<=12)),...
        vphi((vphi<=0 & cbcLT>12)|(vphi>=0 & cbcLT<=12)),'rp')
    scatter(cbcLT((vphiO<=0 & cbcLT>12)|(vphiO>=0 & cbcLT<=12)),...
        vphiO((vphiO<=0 & cbcLT>12)|(vphiO>=0 & cbcLT<=12)),'rp')
    title('Averaged Magnetosheath Flows')
    axis([3 20 -400 400])
end
set(gca,'XTick',[3 6 9 12 15 18])
xlabel('Local Time')
ylabel('V_\phi (km/s)')
set(gca, 'FontSize', 16)

%     figure
%     histogram(grr(grr_LT == 0))
%     title('dusk')
%     xlabel('v_\phi km/s')
%     figure
%     histogram(grr(grr_LT == 1))
%     title('dawn')
%     xlabel('v_\phi km/s')

%     figure
%     hist(crossing_by_crossing_average(cbcLT >= 12 & crossing_by_crossing_average(:,1) ~= 0,1))
%     title('dusk')
%     xlabel('v_\phi km/s')
%     figure
%     hist(crossing_by_crossing_average(cbcLT <= 12 & crossing_by_crossing_average(:,1) ~= 0,1))
%     title('Dawn')
%     xlabel('v_\phi km/s')

    avg_model_data_T = nan(slices,1);
    avg_model_data_density = nan(slices,1);
    avg_model_data_v_phi = nan(slices,1);

    for i = 1:slices
        avg_model_data_density(i) = geomean(crossing_by_crossing_average(floor(k*cbcLT) == i & crossing_by_crossing_average(:,3),3));
        avg_model_data_T(i) =  mean(crossing_by_crossing_average(floor(k*cbcLT) == i & crossing_by_crossing_average(:,4),4));
        avg_model_data_v_phi(i) = mean(crossing_by_crossing_average(floor(k*cbcLT) == i & crossing_by_crossing_average(:,1),1));
    end


    density_first = find(~isnan(avg_model_data_density),1,'first');
    density_last = find(~isnan(avg_model_data_density),1,'last');
    temp_first = find(~isnan(avg_model_data_T),1,'first');
    temp_last = find(~isnan(avg_model_data_T),1,'last');
    v_phi_first = find(~isnan(avg_model_data_v_phi),1,'first');
    v_phi_last = find(~isnan(avg_model_data_v_phi),1,'last');

    interp_density = interp1(find(~isnan(avg_model_data_density)),avg_model_data_density(~isnan(avg_model_data_density)),density_first:density_last);
    interp_temp = interp1(find(~isnan(avg_model_data_T)),avg_model_data_T(~isnan(avg_model_data_T)),temp_first:temp_last);
    interp_vphi = interp1(find(~isnan(avg_model_data_v_phi)),avg_model_data_v_phi(~isnan(avg_model_data_v_phi)),v_phi_first:v_phi_last);

    smooth_density = transpose(smooth(interp_density));
    smooth_temp = transpose(smooth(interp_temp));
    smooth_vphi = transpose(smooth(interp_vphi));

    %smooth_density = interp_density;
    %smooth_temp = interp_temp;
    %smooth_vphi = interp_vphi;

    p1 = polyfit(temp_first:temp_last,smooth_temp,2);
    p2 = polyfit(density_first:density_last,smooth_density,2);
    p22 = polyfit(density_first+4:density_last-6,smooth_density(4:end-7),2);
    p4 = polyfit(v_phi_first:v_phi_last,smooth_vphi,3);
    p44 = polyfit([7 8 9 10 24 38 39 40 41],[-200 -200 -200 -200 0 200 200 200 200],3);

avg_model_data_v_phi(v_phi_first:v_phi_last) = interp_vphi;

    d2 = avg_model_data_v_phi(v_phi_first:22);
    %d22 = avg_model_data_v_phi2(v_phi_first2+4:22);

%LOTS OF PLOTS
%---------------------------------------------------------------------------

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
%--------------------------------------------------------------------------

    temperature = p1;
    density = p22;
    v_r = p3;
    v_phi = p44;
