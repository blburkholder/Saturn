function plot_data(mag_data,sheath_data,SW_data,max_windows)
magsphere_LT_windows = 0;
sheath_inner_LT_windows = 0;
plot_plot_plot = false;
how_many_windows = 20;
h = max_windows/how_many_windows;
time_resolution = 60;
slices = 24*(60/time_resolution);
mult = slices/24;
format shorteng
for w = 1:3
    if w == 1
        f_data = mag_data;
    elseif w == 2
        f_data = sheath_data;
    else
        f_data = SW_data;
    end
    [t,number_of_points] = size(f_data); 

    data = f_data;
    KAW_q_mean = mean(data(1,data(1,:) ~= 0));
    KAW_q_std = std(data(1,data(1,:) ~= 0));
    MHD_q_mean = mean(data(2,data(2,:) ~= 0));
    MHD_q_std = std(data(2,data(2,:) ~= 0));
    KAW_slope_mean = mean(data(9,data(9,:) ~= 0));
    MHD_slope_mean = mean(data(10,data(10,:) ~= 0));
    KAW_slope_std = std(data(9,data(9,:) ~= 0));
    MHD_slope_std = std(data(10,data(10,:) ~= 0));
    K_slope_max = KAW_slope_mean + 2*KAW_slope_std;
    K_slope_min = KAW_slope_mean - 2*KAW_slope_std;
    K_q_max = KAW_q_mean + 2*KAW_q_std;
    M_slope_max = MHD_slope_mean + 2*MHD_slope_std;
    M_slope_min = MHD_slope_mean - 2*MHD_slope_std;
    M_q_max = MHD_q_mean+2*MHD_q_std;

    %filtering
    data(1,data(9,:) > K_slope_max) = 0;
    data(1,data(9,:) < K_slope_min) = 0;
    data(1,data(1,:) > K_q_max) = 0;
    data(2,data(10,:) > M_slope_max) = 0;
    data(2,data(10,:) < M_slope_min) = 0;
    data(2,data(2,:) > M_q_max) = 0;

    color_plot = zeros(180,how_many_windows,2);
    %color_plot_ = ones(180,how_many_windows,2);
    color_plot3 = zeros(slices,how_many_windows,2);
    %color_plot3_ = ones(slices,how_many_windows,2);
    color_plot5 = zeros(75,how_many_windows,2);
    %color_plot5_ = ones(75,how_many_windows,2);
    color_plot7 = zeros(150,how_many_windows,2);
    %color_plot7_ = ones(100,how_many_windows,2);

    %avoid_zeros = (max(f_data(1,:)))^1.5;
    %grossest loop ever
    for y = 1:number_of_points
        LT_index = floor(mult*(f_data(4,y)))+1;

        if f_data(1,y) || f_data(2,y)
            if f_data(1,y)
                M_or_K = 1;
            else
                M_or_K = 2;
            end
            if (90 + f_data(17,y)) > 79 && (90 + f_data(17,y)) < 101
                %latitude
                    color_plot(floor(90 + f_data(17,y)),ceil(f_data(3,y)/h),1) = color_plot(floor(90 + f_data(17,y)),ceil(f_data(3,y)/h),1) + log(f_data(M_or_K,y));
                    color_plot(floor(90 + f_data(17,y)),ceil(f_data(3,y)/h),2) = color_plot(floor(90 + f_data(17,y)),ceil(f_data(3,y)/h),2) + 1;
                %LT
                    color_plot3(LT_index,ceil(f_data(3,y)/h),1) = color_plot3(LT_index,ceil(f_data(3,y)/h),1) + log(f_data(M_or_K,y));
                    color_plot3(LT_index,ceil(f_data(3,y)/h),2) = color_plot3(LT_index,ceil(f_data(3,y)/h),2) + 1;
                %radial
                    color_plot5(floor(f_data(5,y)),ceil(f_data(3,y)/h),1) = color_plot5(floor(f_data(5,y)),ceil(f_data(3,y)/h),1) + log(f_data(M_or_K,y)); 
                    color_plot5(floor(f_data(5,y)),ceil(f_data(3,y)/h),2) = color_plot5(floor(f_data(5,y)),ceil(f_data(3,y)/h),2) + 1;
                %standoff distance
                if f_data(3,y) == 1
                    z = y;
                    dynamic_pressure = get_dynamic_pressure(f_data(5,z),2*pi*(f_data(4,z)/24)-pi);
                    if dynamic_pressure > 0
                        R0 = 10.3*(dynamic_pressure^(-0.2));
                        %in units of half an Rs!!!
                        if R0 - floor(R0) < 0.5
                            index = 2*floor(R0);
                        else
                            index = 2*floor(R0)+1;
                        end
                        if  ~isinf(log(f_data(1,z)))
                            color_plot7(index,1,1) = color_plot7(index,1,1) + log(f_data(1,z));
                            color_plot7(index,1,2) = color_plot7(index,1,2) + 1;
                        elseif ~isinf(log(f_data(2,z)))
                            color_plot7(index,1,1) = color_plot7(index,1,1) + log(f_data(1,2));
                            color_plot7(index,1,2) = color_plot7(index,1,2) + 1;
                        end
                        if (z+1 < number_of_points && f_data(3,z+1) == f_data(3,z) + 1)  || (z > 1 && f_data(3,z-1) == f_data(3,z) + 1)
                            if z+1 < number_of_points && f_data(3,z+1) == f_data(3,z) + 1
                                increment = 1;
                            else
                                increment = - 1;
                            end                    
                            while(z + increment < number_of_points && z + increment > 0 && f_data(3,z+increment) == f_data(3,z) + 1)
                                if ~isinf(log(f_data(1,z + increment)))
                                        color_plot7(index,ceil(f_data(3,z+increment)/h),1) = color_plot7(index,ceil(f_data(3,z+increment)/h),1) + log(f_data(1,z+increment));
                                        color_plot7(index,ceil(f_data(3,z+increment)/h),2) = color_plot7(index,ceil(f_data(3,z+increment)/h),2) + 1;
                                elseif ~isinf(log(f_data(2,z + increment)))
                                        color_plot7(index,ceil(f_data(3,z+increment)/h),1) = color_plot7(index,ceil(f_data(3,z+increment)/h),1) + log(f_data(2,z+increment));
                                        color_plot7(index,ceil(f_data(3,z+increment)/h),2) = color_plot7(index,ceil(f_data(3,z+increment)/h),2) + 1;
                                end
                                z = z + increment;
                            end
                        end
                    end
                end
            end
        end 
    end

    stats1 = color_plot(:,:,2);
    color_plot = color_plot(:,:,1);
    %stats1_ = color_plot_(:,:,2);
    %color_plot_ = color_plot_(:,:,1);
    stats2 = color_plot3(:,:,2);
    color_plot3 = color_plot3(:,:,1);
    %stats2_ = color_plot3_(:,:,2);
    %color_plot3_ = color_plot3_(:,:,1);
    stats3 = color_plot5(:,:,2);
    color_plot5 = color_plot5(:,:,1);
    %stats3_ = color_plot5_(:,:,2);
    %color_plot5_ = color_plot5_(:,:,1);
    stats4 = color_plot7(:,:,2);
    color_plot7 = color_plot7(:,:,1);
    %stats4_ = color_plot7_(:,:,2);
    %color_plot7_ = color_plot7_(:,:,1);

    %compute geometric means
    color_plot(stats1 == 0) = NaN;
    color_plot3(stats2 == 0) = NaN;
    color_plot5(stats3 == 0) = NaN;
    color_plot7(stats4 == 0) = NaN;

    color_plot = exp(color_plot./stats1);
    color_plot3 = exp(color_plot3./stats2);
    color_plot5 = exp(color_plot5./stats3);
    color_plot7 = exp(color_plot7./stats4);

    %normalize!
    color_plot2 = transpose(bsxfun(@rdivide,transpose(color_plot),max(transpose(color_plot)))); 
    color_plot4 = transpose(bsxfun(@rdivide,transpose(color_plot3),max(transpose(color_plot3))));
    color_plot6 = transpose(bsxfun(@rdivide,transpose(color_plot5),max(transpose(color_plot5))));
    color_plot8 = transpose(bsxfun(@rdivide,transpose(color_plot7),max(transpose(color_plot7)))); 

        figure;
        g = pcolor(1:how_many_windows,1:150,log10(color_plot7));
        op = colorbar;
        ylabel(op, 'log(W/m^3)');
        set(g,'edgecolor','none'); 
        title('R_0 = 10.3D_p^{-0.20}');
        grid on

    if plot_plot_plot
        %square plots    
        figure
        pcolor(1:how_many_windows,1:180,log10(color_plot));
        title('Not normalized latitude')

        figure
        pcolor(1:how_many_windows,1:180,log10(color_plot2));
        title('Normalized latitude');

        figure;
        pcolor(1:how_many_windows,1:slices,log10(color_plot3));
        title('not normalized binned by LT');

        figure;
        pcolor(1:how_many_windows,1:slices,log10(color_plot4));
        title('normalized binned by LT');

        figure;
        pcolor(1:how_many_windows,1:75,log10(color_plot5));
        title('not normalized binned by R_s');

        figure;
        pcolor(1:how_many_windows,1:75,log10(color_plot6));
        title('normalized binned by R_s');

        figure;
        pcolor(1:how_many_windows,1:150,log10(color_plot7));
        title('dynamic pressure');

        figure;
        pcolor(1:how_many_windows,1:150,log10(color_plot8));
        title('normalized dynamic pressure');
    end
%polar plots
    if w == 1
        r = (how_many_windows-((1:how_many_windows)-1))'/how_many_windows;
        magsphere_LT_windows = color_plot3;
        region = 'magnetosphere';
        pos = [-1 -1 2 2];
    elseif w == 2
        r = (1:how_many_windows)'/how_many_windows;
        sheath_inner_LT_windows = color_plot3;
        region = 'sheath_{MP}';
        pos = [-1/how_many_windows -1/how_many_windows 2/how_many_windows 2/how_many_windows];
    else
        r = (how_many_windows-((1:how_many_windows)-1))'/how_many_windows;
        sheath_outer_LT_windows = color_plot3;
        region = 'sheath_{BS}';
        pos = [-1 -1 2 2];
    end
    theta = pi+(0:slices-1)*(2*pi/(slices-1));
    t_x = r*cos(theta);
    t_y = r*sin(theta);

    figure
    pcolor(t_x,t_y,log10(transpose(color_plot3)));
    title(region);
    colorbar;
    rectangle('Position',pos,'Curvature',[1 1],'EdgeColor','r','linewidth',3)
    %normalized
    %figure;
    %pcolor(t_x,t_y,log10(transpose(color_plot4)));
    %colorbar;

%line plots
    figure;
    %average q all local times as a function of window number from boundary
    average_q_window = zeros(how_many_windows,1);
    for i = 1:how_many_windows
        average_q_window(i) = log10(geomean(color_plot3(~isnan(color_plot3(:,i)),i)));
    end
    plot(average_q_window)
    title(region)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %standard_error = geostd(log10(line_plot),2)./sqrt(nansum(line_plot,2));
    %std_err_sheath = standard_error;
    %std_err_mag = standard_error;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%average q all window numbers as a function of local time
figure
hold on
combo_plot1 = zeros(slices,1);
combo_plot2 = zeros(slices,1);
combo_plot3 = zeros(slices,1);
for i = 1:slices
    combo_plot1(i) = log10(geomean(magsphere_LT_windows(i,~isnan(magsphere_LT_windows(i,:))),2));
    combo_plot2(i) = log10(geomean(sheath_inner_LT_windows(i,~isnan(sheath_inner_LT_windows(i,:))),2));
    combo_plot3(i) = log10(geomean(sheath_outer_LT_windows(i,~isnan(sheath_outer_LT_windows(i,:))),2));
end
plot(combo_plot1)
plot(fliplr(combo_plot2))
plot(combo_plot3)
%errorbar(geomean(log10(line_plot2),2),std_err_mag);
%errorbar(geomean(log10(line_plot3),2),std_err_sheath);
legend('magnetosphere','sheath inner','sheath outer');
hold off

%FANCIEST PLOT YET
r = ((how_many_windows*3)-((1:(how_many_windows*3))-1))'/(how_many_windows*3);
theta = pi+(0:slices-1)*(2*pi/(slices-1));
t_x = r*cos(theta);
t_y = r*sin(theta);

figure;
g = pcolor(t_x,t_y,log10(transpose(horzcat(sheath_outer_LT_windows,fliplr(sheath_inner_LT_windows),magsphere_LT_windows))));
set(g,'edgecolor','none') 
h = colorbar;
ylabel(h, 'log(W/m^3)')
pos = [-1/3 -1/3 2/3 2/3];
rectangle('Position',pos,'Curvature',[1 1],'EdgeColor','r','linewidth',1)
pos = [-2/3 -2/3 4/3 4/3];
rectangle('Position',pos,'Curvature',[1 1],'EdgeColor','m','linewidth',1)
pos = [-1 -1 2 2];
rectangle('Position',pos,'Curvature',[1 1],'EdgeColor','r','linewidth',1)
grid on
%normalized
%color_plot4 = transpose(bsxfun(@rdivide,transpose(horzcat(fliplr(sheath_LT_windows),...
    %magsphere_LT_windows)),max(transpose(horzcat(fliplr(sheath_LT_windows),magsphere_LT_windows)))));

%figure;
%pcolor(t_x,t_y,log10(transpose(color_plot4)));
%colorbar;
    
%histograms near hot region of sheath (not necessary b/c of error bars?)
%{
    figure
    bar(horzcat(fliplr(sheath_LT_windows(49,:,2)),magsphere_LT_windows(49,:,2)));
    title('12:15');
%}
end