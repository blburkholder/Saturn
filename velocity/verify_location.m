    %data = get_LANL_moments();

    path_Name = '/home/computation/Documents/GitProjects/Saturn/';
    file_name = 'Cassini_orbits_2004-2012.txt';
    formatSpec = ['%d %d %d %d %f %f %f %f %f %f %f %f'];
    data_size = [12 Inf]; 
    full_file_name = horzcat(path_Name, file_name);
    fileID = fopen(full_file_name,'r');
   
    %end = r...... end-1 = LT
    loc_data = fscanf(fileID,formatSpec,data_size);
    fclose(fileID);
    boundaries = get_location_regions_boundary_data();
    [row,col] = size(boundaries);
    boundaries2 = zeros(row-1,col);
    boundaries2(1,:) = boundaries(1,:);
    boundaries2(3:end,:) = boundaries(4:end,:);
    
    for i = 1:col
        d = datevec(strcat(num2str(boundaries(2,i)),'/',num2str(boundaries(3,i)),'/',num2str(boundaries(1,i))));
        v = datenum(d);                
        boundaries2(2,i) = v - datenum(d(1), 1,0);
    end

    location = zeros(5,length(loc_data));
    location(1,:) = datenum(loc_data(1,:),0,loc_data(2,:),loc_data(3,:),0,0);
    location(2:end,:) = loc_data(end-3:end,:);
 
    bounds = zeros(2,col);
    bounds(1,:) = datenum(boundaries2(1,:),0,boundaries2(2,:),boundaries2(3,:),0,0);
    bounds(2,:) = boundaries2(end,:);

figure
hold on
    for i = 1:col-1
            trav = location(2:end,location(1,:) >= bounds(1,i) & location(1,:) <= bounds(1,i+1));
            r = trav(4,:);
            theta = pi*(trav(3,:))/12 + pi;
            x = r.*cos(theta);
            y = r.*sin(theta);
            if bounds(2,i) == 2
                plot(x,y,'r');
            elseif bounds(2,i) == 1 || bounds(2,i) == -2
                plot(x,y,'g');
            else
                plot(x,y,'b');
            end
        end
   

    %x = loc_data(end,:).*cos(loc_data(end-1,:));
    %y = loc_data(end,:).*sin(loc_data(end-1,:));
    %plot(x,y)
    
    %{
    count  = 0;
    for i = 1:length(data)
        a = loc_data(end-3:end,loc_data(1,:) == data(2,i) & loc_data(2,:) == data(3,i) & loc_data(3,:) == data(4,i));
        if ~isempty(a)
            r = sqrt((a(4) - data(26,i))^2);
            lat = sqrt((a(1) - data(27,i))^2);
            lt = sqrt((a(3) - data(28,i))^2);
            if (lt > 1 & lt < 23) | r > 1 | lat > 5
                data(2,i)
                data(3,i)
                data(4,i)
                pause
            end
        end
    end
%}