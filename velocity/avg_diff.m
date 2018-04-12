avg_dif = zeros(10,1);
avg_dif(1) = abs(abs(means_time3(1,2)) - abs(means_time3(1,1)));
for i = 2:10
    avg_dif(i) = ((i-1)*avg_dif(i-1) + abs(abs(means_time3(i,2)) - abs(means_time3(i,1))))/i;
end
figure
plot(avg_dif)