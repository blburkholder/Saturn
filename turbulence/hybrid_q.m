%x = ans.X;
bx = ans.BX;
by = ans.BY;
bz = ans.BZ;
ux = ans.UX;
uz = ans.UZ;
rhon = mean(mean(ans.DENS));
%temp = ans.TEMP_P;

bx = bx(:,2:end-1);
by = by(:,2:end-1);
bz = bz(:,2:end-1);

R = 25;
lat = 0;
avagdro_number = 6.022140857e23; 
ion_mass = 1e-3/avagdro_number;
mu_0 = 4*pi()*1e-7; %N/A^2

%[v_phi,v_r] = get_v_rel(R);
v_phi = mean(mean(ux));
v_r = mean(mean(uz));
H = get_scale_height(R);
%rhon = get_density(lat,H,R)
rho = rhon*ion_mass;
temp = get_temperature(H);

B_vector_mean = get_B_vector_mean(mean(mean(bx)), mean(mean(by)), mean(mean(bz)));

[ll,height] = size(bx);
B_vector = zeros(3,ll,height);
B_mag = zeros(3,ll,height);
B_vector(1,:,:) = bx;
B_vector(2,:,:) = by;
B_vector(3,:,:) = bz;
    
B_ave(1,1:ll,1:height) = B_vector_mean(1);
B_ave(2,1:ll,1:height) = B_vector_mean(2);
B_ave(3,1:ll,1:height) = B_vector_mean(3);
    
unit_B_ave = zeros(3,ll,height);
unit_B_ave(1,:,:) = B_ave(1,:,:)./sqrt(B_ave(1,:,:).^2 + B_ave(2,:,:).^2 + B_ave(3,:,:).^2);
unit_B_ave(2,:,:) = B_ave(2,:,:)./sqrt(B_ave(1,:,:).^2 + B_ave(2,:,:).^2 + B_ave(3,:,:).^2);
unit_B_ave(3,:,:) = B_ave(3,:,:)./sqrt(B_ave(1,:,:).^2 + B_ave(2,:,:).^2 + B_ave(3,:,:).^2);
    
B_pertr = B_vector - B_ave;
    
B_parallel_pertr(1,:,:) = dot(B_pertr, unit_B_ave).*unit_B_ave(1,:,:);
B_parallel_pertr(2,:,:) = dot(B_pertr, unit_B_ave).*unit_B_ave(2,:,:);
B_parallel_pertr(3,:,:) = dot(B_pertr, unit_B_ave).*unit_B_ave(3,:,:);
    
B_perp_pertr = B_pertr - B_parallel_pertr;

figure; hold on
bhist = sum(B_perp_pertr.^2);
bbhist = sum(B_parallel_pertr.^2);
bbbhist = sum(B_pertr.^2);
histogram(log10(bhist(:)))
histogram(log10(bbhist(:)))
histogram(log10(bbbhist(:)))

lambda = 360000;                 %grid separation (in ion inertial lengths?)
ks = 2*pi/lambda;            % sampling wavenumber                      
L1 = 129; % Length of signal (in ion intertial lengths?)
%L2 = 137;
L2 = 108;

perpX = fft2(squeeze(B_perp_pertr(1,:,:)));
perpY = fft2(squeeze(B_perp_pertr(2,:,:)));
perpZ = fft2(squeeze(B_perp_pertr(3,:,:)));
[s1,s2] = size(perpZ);
k1 = ks*(-L1/2+1:L1/2)/L1;
k2 = ks*(-L2/2+1:L2/2)/L2;
[kgrid1,kgrid2] = meshgrid(k2,k1);

pii = (0:(pi/65):(pi))/lambda;
sumsum = 0;
Bpar2_B02 = 1e-18/(B_vector_mean(1)^2+B_vector_mean(2)^2+B_vector_mean(3)^2);

%figure
for i = 1:65
cond1 = ((kgrid1-k1(ceil(L1/2))).^2+(kgrid2-k2(ceil(L2/2))).^2) > pii(i)^2  & ((kgrid1-k1(ceil(L1/2))).^2+(kgrid2-k2(ceil(L2/2))).^2) < pii(i+1)^2;
absx = abs(fftshift(perpX));
absz = abs(fftshift(perpZ));
absy = abs(fftshift(perpY));
fluct = ((absx.^2+absy.^2+absz.^2)*pii(i)).^(3/2);
sumsum = sumsum + pii(i)*sum(sum(fluct(cond1)))/sqrt(rho*mu_0^3);

% absx(cond1) = NaN;
% absz(cond1) = NaN;
% pcolor(k2,k1,log10(absx.^2+absz.^2))
% shading interp
% colorbar
% xlabel('k_{\perpx}')
% ylabel('k_{\perpy}')
%     colorbar
%     pause(0.1)
end

q = sumsum/129

bsqr = (B_vector_mean(1)^2+B_vector_mean(2)^2+B_vector_mean(3)^2);
mass_density = rho;
beta = 1;
kz = 2*pi/(60000000);
emoob = exp(-1/beta);
Va = sqrt(bsqr)/sqrt(mu_0*mass_density);
d_perp = sqrt(pi/2)*Bpar2_B02*Va/(kz)*emoob
