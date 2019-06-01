% load(['q_d_data/lat.mat'])
% load(['q_d_data/LT.mat'])
% load(['q_d_data/R'])
% load(['q_d_data/tfb'])
% load(['q_d_data/fucktuation_perp'])
% load(['q_d_data/fucktuation_par'])
% load(['q_d_data/mva_eigvals.mat'])
% load(['q_d_data/nl_measure.mat'])
% load(['q_d_data/q'])
% load(['q_d_data/d_complicated_integral'])
% load(['q_d_data/d_beta_exponential'])
% load(['q_d_data/current_sheet'])

load(['lat.mat'])
load(['LT.mat'])
load(['R'])
load(['tfb'])
load(['fucktuation_perp'])
load(['fucktuation_par'])
load(['mva_eigvals.mat'])
load(['nl_measure.mat'])
load(['q'])
load(['d_complicated_integral'])
load(['d_beta_exponential'])
load(['current_sheet'])

d = d_beta_exponential;

fucktuation = sqrt(fucktuation_perp.^2+fucktuation_par.^2);
%q = q/sqrt(0.001);
%d = d/sqrt(0.001);

figure
scatter(LT(~isnan(q)),log10(mva_eigvals(~isnan(q))),[],log10(q(~isnan(q))),'.');
ylabel('$\mathcal{E}$ (T)','interpreter','latex')
set(gca, 'FontSize', 16)
xlabel('Local Time')
title('variance')
h = colorbar;
set(get(h,'title'),'string','q (W/m^3)')
colormap(jet)
set(gca,'XTick',[3 6 9 12 15 18])
axis([3 20 -23 -17])

figure
scatter(LT(~isnan(q)),log10(nl_measure(~isnan(q))),[],log10(q(~isnan(q))),'.');
title('embedding')
set(gca, 'FontSize', 16)
xlabel('Local Time')
ylabel('$\mathcal{N}$','interpreter','latex')
h = colorbar;
set(get(h,'title'),'string','q (W/m^3)')
colormap(jet)
set(gca,'XTick',[3 6 9 12 15 18])
axis([3 20 -7.5 -5])

figure
scatter(LT(~isnan(q)),log10(fucktuation(~isnan(q))),[],log10(q(~isnan(q))),'.');
title('fluctuation')
set(gca, 'FontSize', 16)
xlabel('Local Time')
ylabel('\delta B (T)')
h = colorbar;
set(get(h,'title'),'string','q (W/m^3)')
colormap(jet)
set(gca,'XTick',[3 6 9 12 15 18])
axis([3 20 -11 -8])

figure
scatter(LT(~isnan(d)),log10(mva_eigvals(~isnan(d))),[],log10(d(~isnan(d))),'.');
title('variance')
set(gca, 'FontSize', 16)
xlabel('Local Time')
ylabel('$\mathcal{E}$ (T)','interpreter','latex')
h = colorbar;
h.Label.String = 'D_\perp (m^2/s)';
set(get(h,'title'),'string','D_\perp (m^2/s)')
colormap(jet)
set(gca,'XTick',[3 6 9 12 15 18])
axis([3 20 -23 -17])
%caxis([2 log10(max(d))])

figure
scatter(LT(~isnan(d)),log10(nl_measure(~isnan(d))),[],log10(d(~isnan(d))),'.');
title('embedding')
set(gca, 'FontSize', 16)
xlabel('Local Time')
ylabel('$\mathcal{N}$','interpreter','latex')
h = colorbar;
set(get(h,'title'),'string','D_\perp (m^2/s)')
colormap(jet)
set(gca,'XTick',[3 6 9 12 15 18])
axis([3 20 -7.5 -5])
%caxis([2 log10(max(d))])

figure
scatter(LT(~isnan(d)),log10(fucktuation(~isnan(d))),[],log10(d(~isnan(d))),'.');
title('fluctuation')
set(gca, 'FontSize', 16)
xlabel('Local Time')
ylabel('\delta B (T)')
h = colorbar;
set(get(h,'title'),'string','D_\perp (m^2/s)')
colormap(jet)
set(gca,'XTick',[3 6 9 12 15 18])
axis([3 20 -11 -8])
%caxis([2 log10(max(d))])

figure
%h = histogram(log10(nl_measure),'Normalization','probability');
h = histogram(log10(nl_measure));
heights = h.Values;
xv = h.BinEdges;
xvals = ( xv(2:end)+xv(1:end-1) )/2;
hold on
f = fit(xvals',heights','gauss2');
plot(f,xvals,heights)
xx = xvals(1):0.01:xvals(end);
yvals = f.a1*exp(-((xx-f.b1)/f.c1).^2) + f.a2*exp(-((xx-f.b2)/f.c2).^2);
ind = find(yvals == max(yvals));
stdd = max([f.c1,f.c2]);
plot([xx(ind)+stdd,xx(ind)+stdd],[0 max(yvals)])
thresh1 = xx(ind) + stdd;
title(['nonlinear measure - ',num2str(100*sum(log10(nl_measure) >= thresh1 & current_sheet)/length(nl_measure)),'% active'])
set(gca, 'FontSize', 12)

figure
h = histogram(log10(mva_eigvals));
heights = h.Values;
xv = h.BinEdges;
xvals = ( xv(2:end)+xv(1:end-1) )/2;
hold on
f = fit(xvals',heights','gauss2');
plot(f,xvals,heights)
xx = xvals(1):0.01:xvals(end);
yvals = f.a1*exp(-((xx-f.b1)/f.c1).^2) + f.a2*exp(-((xx-f.b2)/f.c2).^2);
ind = find(yvals == max(yvals));
stdd = max([f.c1,f.c2]);
plot([xx(ind)+stdd,xx(ind)+stdd],[0 max(yvals)])
thresh2 = xx(ind) + stdd;
title(['variance - ',num2str(100*sum(log10(mva_eigvals) >= thresh2 & current_sheet)/length(mva_eigvals)),'% active'])
set(gca, 'FontSize', 12)

figure
h = histogram(log10(fucktuation));
heights = h.Values;
xv = h.BinEdges;
xvals = ( xv(2:end)+xv(1:end-1) )/2;
hold on
f = fit(xvals',heights','gauss2');
plot(f,xvals,heights)
xx = xvals(1):0.01:xvals(end);
yvals = f.a1*exp(-((xx-f.b1)/f.c1).^2) + f.a2*exp(-((xx-f.b2)/f.c2).^2);
ind = find(yvals == max(yvals));
stdd = max([f.c1,f.c2]);
plot([xx(ind)+stdd,xx(ind)+stdd],[0 max(yvals)])
thresh3 = xx(ind) + stdd;
title(['fluctuation - ',num2str(100*sum(log10(fucktuation) >= thresh3 & current_sheet)/length(fucktuation)),'% active'])
set(gca, 'FontSize', 12)

figure
x = R.*cos(pi*LT/12-pi);
y = R.*sin(pi*LT/12-pi);
scatter(x,y,'k.')
hold on
plot([-60 40],[0 0],'k')
scatter(x(log10(nl_measure) >= thresh1 & current_sheet),y(log10(nl_measure) >= thresh1 & current_sheet),[150],'c.')
scatter(x(log10(mva_eigvals) >= thresh2 & current_sheet),y(log10(mva_eigvals) >= thresh2 & current_sheet),[150],'m.')
scatter(x(log10(fucktuation) >= thresh3 & current_sheet),y(log10(fucktuation) >= thresh3 & current_sheet),[150],'y.')
%scatter(x(log10(nl_measure) >= thresh1 & log10(mva_eigvals) >= thresh2 & log10(fucktuation) >= thresh3),...
    %y(log10(nl_measure) >= thresh1 & log10(mva_eigvals) >= thresh2 & log10(fucktuation) >= thresh3) ,'rp')
title(['active ',num2str(100*sum((log10(nl_measure) >=...
    thresh1 | log10(mva_eigvals) >= thresh2 | log10(fucktuation) >=...
    thresh3) & current_sheet)/length(R)),'% \pm',num2str(100/sqrt(4*length(R))),'%'])
legend('all windows','reference line','embedding','variance','fluctuation','AND')
ylabel('Y (R_s)')
xlabel('X (R_s)')

qtime = zeros(1,10);
dtime = zeros(1,10); 
d1time = zeros(1,10);
qtime_err = zeros(1,10);
dtime_err = zeros(1,10);
d1time_err = zeros(1,10);
for i = 1:10
    qtime(i) = geomean(q(tfb > (i-1)*5 & tfb < i*5 & ~isnan(q)));
    dtime(i) = geomean(d_complicated_integral(tfb > (i-1)*5 & tfb < i*5 & ~isnan(d_complicated_integral) & d_complicated_integral ~= 0));
    d1time(i) = geomean(d(tfb > (i-1)*5 & tfb < i*5 & ~isnan(d)));
    qtime_err(i) = std(log10(q(tfb > (i-1)*5 & tfb < i*5 & ~isnan(q))))/sqrt(sum(tfb > (i-1)*5 & tfb < i*5 & ~isnan(q)));
    dtime_err(i) = std(log10(d_complicated_integral(tfb > (i-1)*5 & tfb < i*5 & ~isnan(d_complicated_integral) & d_complicated_integral ~= 0)))/sqrt(sum(tfb > (i-1)*5 & tfb < i*5 & ~isnan(d_complicated_integral& d_complicated_integral ~= 0)));
    d1time_err(i) = std(log10(d(tfb > (i-1)*5 & tfb < i*5 & ~isnan(d))))/sqrt(sum(tfb > (i-1)*5 & tfb < i*5 & ~isnan(d)));
end

figure
histogram(log10(q(LT>12 & ~isnan(q))),'facecolor','green','BinWidth',0.1)
hold on
histogram(log10(q(LT<=12 & ~isnan(q))),'facecolor','blue','BinWidth',0.1)
xlabel('log_{10}(q)')
set(gca, 'FontSize', 24)
histogram(log10(q((log10(nl_measure)>=thresh1 | log10(mva_eigvals)>=thresh2 | log10(fucktuation) >= thresh3) & current_sheet)),'facecolor','red','BinWidth',0.1)
box off;
axis([-21 -12 0 1400])
legend('dusk','dawn','active')
ax1 = gca; % current axes
ax1.XColor = 'k';
ax1.YColor = 'k';
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none','XColor','r','YColor','r');
line(10:10:100,log10(qtime),'Parent',ax2,'Color','r','LineWidth',2)
hold on
errorbar(10:10:100,log10(qtime),qtime_err,'Parent',ax2,'Color','r','LineWidth',2)
xlabel('time from boundary')
ylabel('log_{10}(q)')
set(gca, 'FontSize', 24)

figure
histogram(log10(d_complicated_integral(LT >12)),'facecolor','green','BinWidth',0.1)
hold on
histogram(log10(d_complicated_integral(LT <= 12)),'facecolor','blue','BinWidth',0.1)
xlabel('log_{10}(D_\perp)')
set(gca, 'FontSize', 24)
histogram(log10(d_complicated_integral((log10(nl_measure)>=thresh1 | log10(mva_eigvals)>=thresh2 | log10(fucktuation) >= thresh3) & current_sheet)),'facecolor','red','BinWidth',0.1)
box off;
%axis([0 11 0 1600])
legend('dusk','dawn','active')
ax1 = gca; % current axes
ax1.XColor = 'k';
ax1.YColor = 'k';
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none','XColor','r','YColor','r');
line(10:10:100,log10(dtime),'Parent',ax2,'Color','r','LineWidth',2)
hold on
errorbar(10:10:100,log10(dtime),dtime_err,'Parent',ax2,'Color','r','LineWidth',2)
xlabel('time from boundary')
ylabel('log_{10}(D_\perp)')
set(gca, 'FontSize', 24)

figure
histogram(log10(d_beta_exponential(LT >12)),'facecolor','green','BinWidth',0.1)
hold on
histogram(log10(d_beta_exponential(LT <= 12)),'facecolor','blue','BinWidth',0.1)
xlabel('log_{10}(D_\perp)')
set(gca, 'FontSize', 24)
histogram(log10(d_beta_exponential((log10(nl_measure)>=thresh1 | log10(mva_eigvals)>=thresh2 | log10(fucktuation) >= thresh3) & current_sheet)),'facecolor','red','BinWidth',0.1)
box off;
%axis([0 11 0 1600])
legend('dusk','dawn','active')
ax1 = gca; % current axes
ax1.XColor = 'k';
ax1.YColor = 'k';
ax1_pos = ax1.Position; % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none','XColor','r','YColor','r');
line(10:10:100,log10(d1time),'Parent',ax2,'Color','r','LineWidth',2)
hold on
errorbar(10:10:100,log10(d1time),d1time_err,'Parent',ax2,'Color','r','LineWidth',2)
xlabel('time from boundary')
ylabel('log_{10}(D_\perp)')
set(gca, 'FontSize', 24)