function [] = myflag_vs_michelleflag(flagged_file)

    data = get_LANL_moments();
    dates = 24*60*(datenum(data(2,:),1,1) + (data(3,:)-1) + data(4,:)/24 + data(5,:)/(24*60)...
            + data(6,:)/(24*60*60) - datenum(2004,1,1));

    fid = fopen(strcat('/home/computation/Documents/GitProjects/',flagged_file));
    FC = textscan(fid,'%s','Delimiter','\n');
    a = size(FC{1});
    flag_view = zeros(1,a(1)-1);
    flag_edge_anode = zeros(1,a(1)-1);
    flagged_date = zeros(1,a(1)-1);
    LT = zeros(1,a(1)-1);
    for i = 2:a(1)
        b = FC{1}(i);
        flag_view(i-1) = b{1}(1);
        flag_edge_anode(i-1) = b{1}(2);
        date_split = strsplit(b{1});
        year = date_split(1);
        flagged_date(i-1) = 24*60*(datenum(str2double((year{1}(3:end))),1,1) + (str2double(cell2mat(date_split(2)))-1) +...
            str2double(cell2mat(date_split(3)))/24 + str2double(cell2mat(date_split(4)))/(24*60)...
            + str2double(cell2mat(date_split(5)))/(24*60*60) - datenum(2004,1,1));
        LT(i-1) = str2double(cell2mat(date_split(6)));
    end
   
    vphi_data = zeros(1,length(flagged_date));
    vphi_data_flag = zeros(1,length(flagged_date));
    for k =1:length(flagged_date)
        vphi_data(k) = data(10,dates == flagged_date(k));
        vphi_data_flag(k) = data(37,dates == flagged_date(k));
    end

    Michelle_vphi = vphi_data((vphi_data_flag & LT < 12) | (~vphi_data_flag & LT >= 12));
    Michelle_LT = LT((vphi_data_flag & LT < 12) | (~vphi_data_flag & LT >= 12));
    My_vphi = vphi_data(flag_view == 105 & flag_edge_anode == 43);
    My_LT  = LT(flag_view == 105 & flag_edge_anode == 43);

    figure
    hold on
    scatter([1,2,3,4,5,6,7,8],[mean(Michelle_vphi(Michelle_LT < 10.5)),mean(Michelle_vphi(Michelle_LT < 11 & Michelle_LT >= 10.5)),...
        mean(Michelle_vphi(Michelle_LT < 11.5 & Michelle_LT >= 11)),mean(Michelle_vphi(Michelle_LT < 12 & Michelle_LT >= 11.5)),...
        mean(Michelle_vphi(Michelle_LT < 12.5 & Michelle_LT >= 12)),mean(Michelle_vphi(Michelle_LT < 13 & Michelle_LT >= 12.5)),...
        mean(Michelle_vphi(Michelle_LT < 13.5 & Michelle_LT >= 13)),mean(Michelle_vphi(Michelle_LT <= 14 & Michelle_LT >= 13.5))])
    scatter([1,2,3,4,5,6,7,8],[mean(My_vphi(My_LT < 10.5)),mean(My_vphi(My_LT < 11 & My_LT >= 10.5)),...
        mean(My_vphi(My_LT < 11.5 & My_LT >= 11)),mean(My_vphi(My_LT < 12 & My_LT >= 11.5)),...
        mean(My_vphi(My_LT < 12.5 & My_LT >= 12)),mean(My_vphi(My_LT < 13 & My_LT >= 12.5)),...
        mean(My_vphi(My_LT < 13.5 & My_LT >= 13)),mean(My_vphi(My_LT <= 14 & My_LT >= 13.5))],'*')


end