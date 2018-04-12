h =  findobj('type','figure');
%t = [10,11,19,20];
for i = 1:length(h)
    saveas(h(i),strcat('assorted_fig',num2str(i),'.png'))
    %saveas(h(i),strcat('assorted_fig',num2str(t(i)),'.png'))
end