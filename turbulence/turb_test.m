notzero = allmodel_10mins_msphere(1,:) ~= 0;
data = allmodel_10mins_msphere(:,notzero);

winds = 10;
TF = 24;

turb = zeros(winds,TF);

for i  = 1:winds
    these = data(:,data(3,:) == i);
    if sum(these(3,:) == i) > 0
        for j =1:TF
            those = these(:,ceil(these(4,:)) == j);
            turb(i,j) = geomean(those(1,:));
        end
    end
end
% 
% r = winds:-1:1;
% theta = pi*(1:TF)/12-pi;
% 
% [r,theta] = meshgrid(r,theta);
% 
% figure
% pcolor(r'.*cos(theta'),r'.*sin(theta'),log10(turb))
% colorbar

sizes = -log10(data(1,:)) + max(log10(data(1,:)));
sizes = max(sizes) - sizes+0.01;

figure
scatter(data(4,:),data(3,:)/2,500*sizes,log10(data(1,:)),'.')
hold on
scatter(0,1)
scatter(0,1)
scatter(0,1)
scatter(0,1)
scatter(0,1)
scatter(0,1)
scatter(0,1)
scatter(0,1)
scatter(0,1)

colorbar
tmean = nanmean(turb');
legend(num2str(tmean(1)),num2str(tmean(2)),...
    num2str(tmean(3)),num2str(tmean(4)),num2str(tmean(5)),...
    num2str(tmean(6)),num2str(tmean(7)),num2str(tmean(8)),...
    num2str(tmean(9)),num2str(tmean(10)))

figure
hold on
for i = 1:10
    q = data(1,ceil(data(3,:)) == i);
    LT = data(4,ceil(data(3,:)) == i);
    %wind = data(3,ceil(data(3,:)) == i);
    scatter(LT,log10(q),'.')
end



