[location_data] = get_location_data(); 
rbd = get_location_regions_boundary_data();

diffs = datenum(rbd(1,2:end),rbd(2,2:end),rbd(3,2:end),rbd(4,2:end),rbd(5,2:end),rbd(6,2:end))...
        - datenum(rbd(1,1:end-1),rbd(2,1:end-1),rbd(3,1:end-1),rbd(4,1:end-1),rbd(5,1:end-1),rbd(6,1:end-1));
rx = zeros(length(diffs)+1,1);

for i = 1:length(diffs)
    rx(i+1) = rx(i) + diffs(i);
end

active_date1 = window_dates(window_dates>0)/(24*60)+datenum(2004,1,1)-...
    datenum(rbd(1,1),rbd(2,1),rbd(3,1),rbd(4,1),rbd(5,1),rbd(6,1));
% active_date1 = nl_measure_active_dates(nl_measure_active_dates>0)-...
%     datenum(rbd(1,1),rbd(2,1),rbd(3,1),rbd(4,1),rbd(5,1),rbd(6,1));
% active_date2 = mva_eigvals_active_dates(mva_eigvals_active_dates>0)-...
%     datenum(rbd(1,1),rbd(2,1),rbd(3,1),rbd(4,1),rbd(5,1),rbd(6,1));
%active_date = starting_active_time/(24*60) -...
%    datenum(rbd(1,1),rbd(2,1),rbd(3,1),rbd(4,1),rbd(5,1),rbd(6,1));

figure
hold on
subplot(8,1,1)
for i = 1:length(rx)-1
    time_line(rbd,rx,i)
    if rx(i+1) > 2625
        subplot(8,1,8)
        time_line(rbd,rx,i)
    elseif rx(i+1) > 2250
        subplot(8,1,7)
        time_line(rbd,rx,i)
    elseif rx(i+1) > 1875
        subplot(8,1,6)
        time_line(rbd,rx,i)
    elseif  rx(i+1) > 1500
        subplot(8,1,5)
        time_line(rbd,rx,i)
    elseif rx(i+1) > 1125
        subplot(8,1,4)
        time_line(rbd,rx,i)
    elseif  rx(i+1) > 750
        subplot(8,1,3)
        time_line(rbd,rx,i)
    elseif rx(i+1) > 375
        subplot(8,1,2)
        time_line(rbd,rx,i)
    end
end

s = 200;
for i = 1:length(active_date1)-1
        if active_date1(i) < 375
            subplot(8,1,1)
        elseif active_date1(i) < 750
            subplot(8,1,2)
        elseif active_date1(i) < 1125
            subplot(8,1,3)
        elseif active_date1(i) < 1500
            subplot(8,1,4)
        elseif active_date1(i) < 1875
            subplot(8,1,5)
        elseif active_date1(i) < 2260
            subplot(8,1,6)
        elseif active_date1(i) < 2625
            subplot(8,1,7)
        else
            subplot(8,1,8)
        end
        h = patch([active_date1(i)-s/(24*60),active_date1(i)-s/(24*60),...
        active_date1(i)+s/(24*60),active_date1(i)+s/(24*60)],[1.0,2.0,2.0,1.0],'m');
        set(h,'EdgeColor','none');
end
        h = patch([active_date1(end)-s/(24*60),active_date1(end)-s/(24*60),...
        active_date1(end)+s/(24*60),active_date1(end)+s/(24*60)],[1.0,2.0,2.0,1.0],'m');
        set(h,'EdgeColor','none');

% for i = 2:length(active_date2)-1
%         if active_date2(i) < 375 
%             subplot(8,1,1)
%         elseif active_date2(i) < 750 
%             subplot(8,1,2)
%         elseif active_date2(i) < 1125 
%             subplot(8,1,3)
%         elseif active_date2(i) < 1500 
%             subplot(8,1,4)
%         elseif active_date2(i) < 1875 
%             subplot(8,1,5)
%         elseif active_date2(i) < 2250 
%             subplot(8,1,6)
%         elseif active_date2(i) < 2625
%             subplot(8,1,7)
%         else
%             subplot(8,1,8)
%         end
%         h = patch([active_date2(i)-s/(24*60),active_date2(i)-s/(24*60),...
%         active_date2(i)+s/(24*60),active_date2(i)+s/(24*60)],[3.0,2.0,2.0,3.0],'c');
%         set(h,'EdgeColor','none');
% end
%         h = patch([active_date2(end)-s/(24*60),active_date2(end)-s/(24*60),...
%         active_date2(end)+s/(24*60),active_date2(end)+s/(24*60)],[3.0,2.0,2.0,3.0],'c');
%         set(h,'EdgeColor','none');

subplot(8,1,1)
axis([0 375 -1 3])
subplot(8,1,2)
axis([375 750 -1 3])
subplot(8,1,3)
axis([750 1125 -1 3])
subplot(8,1,4)
axis([1125 1500 -1 3])
subplot(8,1,5)
axis([1500 1875 -1 3])
subplot(8,1,6)
axis([1875 2250 -1 3])
subplot(8,1,7)
axis([2250 2625 -1 3])
subplot(8,1,8)
axis([2625 3000 -1 3])
xlabel('Cassini days since arrival at Saturn')

% ax = gca;
% axis([rx(1) rx(end) -2.5 2.5])
% ax.YTick = [-2 -1 0 1 2];
% ax.YTickLabel = [{'plasma sheet','solar wind','magnetosheath','magnetosphere','VA active'}];
% xlabel('Cassini days since start of mission');
% daspect([100 1 1]);
