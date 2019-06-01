x = ans.X;
bx = ans.BX;
by = ans.BY;
bz = ans.BZ;
ux = ans.UX;
uz = ans.UZ;
rhon = mean(mean(ans.DENS));
%temp = ans.TEMP_P; %im not sure what units this is so im not using it
%because it only goes into calculating rho_i

%the edges are no fun
bx = bx(:,2:end-1);
by = by(:,2:end-1);
bz = bz(:,2:end-1);

%some of this not used
R = 25;
lat = 0;
avagdro_number = 6.022140857e23; 
ion_mass = 1e-3/avagdro_number;
mu_0 = 4*pi()*1e-7; %N/A^2
kB = 1.3806488e-23; %J/K

%getting parameters for calculation
%[v_phi,v_r] = get_v_rel(R);
v_phi = mean(mean(ux));
v_r = mean(mean(uz));
H = get_scale_height(R);
%rhon = get_density(lat,H,R)
rho = rhon*ion_mass;
%use this because it is the correct units for calculating rho_i
temp = get_temperature(H);

B_vector_mean = get_B_vector_mean(mean(mean(bx)), mean(mean(by)), mean(mean(bz)));
bsqr = (B_vector_mean(1)^2+B_vector_mean(2)^2+B_vector_mean(3)^2);

[ll,height] = size(bx);
B_vector = zeros(3,ll,height);
B_mag = zeros(3,ll,height);
B_parallel_pertr = zeros(3,ll,height);
B_ave = zeros(3,ll,height);
B_vector(1,:,:) = bx;
B_vector(2,:,:) = by;
B_vector(3,:,:) = bz;

%calculate parallel and perp fluctuation
B_ave(1,:,:) = B_vector_mean(1);
B_ave(2,:,:) = B_vector_mean(2);
B_ave(3,:,:) = B_vector_mean(3);
    
unit_B_ave = zeros(3,ll,height);
unit_B_ave(1,:,:) = B_ave(1,:,:)./sqrt(B_ave(1,:,:).^2 + B_ave(2,:,:).^2 + B_ave(3,:,:).^2);
unit_B_ave(2,:,:) = B_ave(2,:,:)./sqrt(B_ave(1,:,:).^2 + B_ave(2,:,:).^2 + B_ave(3,:,:).^2);
unit_B_ave(3,:,:) = B_ave(3,:,:)./sqrt(B_ave(1,:,:).^2 + B_ave(2,:,:).^2 + B_ave(3,:,:).^2);
    
B_pertr = B_vector - B_ave;
    
B_parallel_pertr(1,:,:) = dot(B_pertr, unit_B_ave).*unit_B_ave(1,:,:);
B_parallel_pertr(2,:,:) = dot(B_pertr, unit_B_ave).*unit_B_ave(2,:,:);
B_parallel_pertr(3,:,:) = dot(B_pertr, unit_B_ave).*unit_B_ave(3,:,:);
    
B_perp_pertr = B_pertr - B_parallel_pertr;

% figure; hold on
% bhist = sum(B_perp_pertr.^2);
% bbhist = sum(B_parallel_pertr.^2);
% bbbhist = sum(B_pertr.^2);
% histogram(log10(bhist(:)))
% histogram(log10(bbhist(:)))
% histogram(log10(bbbhist(:)))

lambda = x(2)-x(1);                 %grid separation (in ion inertial lengths?) %need to change me for different resolutions
ks = 2*pi/lambda;            % sampling wavenumber 
   
[L1,L2] = size(bx);
                  
perpX = fft2(squeeze(B_perp_pertr(1,:,:))); 
perpY = fft2(squeeze(B_perp_pertr(2,:,:))); 
perpZ = fft2(squeeze(B_perp_pertr(3,:,:))); 
parX = fft2(squeeze(B_parallel_pertr(1,:,:))./(sqrt(bsqr)));  
parY = fft2(squeeze(B_parallel_pertr(2,:,:))./sqrt((bsqr))); 
parZ = fft2(squeeze(B_parallel_pertr(3,:,:))./sqrt((bsqr)));
[s1,s2] = size(perpZ);
k1 = ks*(0:(1/L1):1/2);
k2 = ks*(0:(1/L2):1/2);
[kgrid1,kgrid2] = meshgrid([-fliplr(k1),k1(2:end)],[-fliplr(k2),k2(2:end)]);

%power spectral density
absx = fftshift(perpX); abs1x = absx.*conj(absx);
absy = fftshift(perpZ); abs1y = absy.*conj(absy);
absz = fftshift(perpY); abs1z = absz.*conj(absz);
absx = fftshift(parX); abs2x = absx.*conj(absx);
absy = fftshift(parZ); abs2y = absy.*conj(absy);
absz = fftshift(parY); abs2z = absz.*conj(absz);
%P(k_perp)

L = L1;
k = k1;
for i = 1:floor(L/2)
cond1 = (kgrid1.^2+kgrid2.^2) >= k(i)^2  & (kgrid1.^2+kgrid2.^2) <= k(i+1)^2;

pii(i) = (k(i)+k(i+1))/2; 

power_spectrum_perp = mean(mean(sqrt(abs1x(cond1').^2+abs1y(cond1').^2+abs1z(cond1').^2)));
fluct_perp = (pii(i)*power_spectrum_perp)^(3/2);
sumsum(i) = pii(i)*fluct_perp/sqrt(rho*mu_0^3);

power_spectrum_par = mean(mean(sqrt(abs2x(cond1').^2+abs2y(cond1').^2+abs2z(cond1').^2)));
fluct_par = power_spectrum_par*pii(i);
Bpar2_B02(i) = fluct_par;

% abs1x(cond1') = NaN;
% abs1z(cond1') = NaN;
% pcolor(log10(abs1x.^2+abs1z.^2)')
% shading interp
% colorbar
% xlabel('k_{\perpx}')
% ylabel('k_{\perpy}')
% daspect([1 1 1])
%     colorbar
%     pause(0.1)
end

larmori = get_larmor_radius(temp,sqrt(bsqr),0.018);
q = nanmean(sumsum)

Ti = 500;
Te = 200;
Va = sqrt(bsqr)/sqrt(mu_0*rho);
%[j1,j4] = evalintegral(sqrt(bsqr),Ti,Te,pii.*larmori,Bpar2_B02,larmori);
kz = 2*pi/(60000000);
krho = pii*larmori;
d_perp = sqrt(pi/8)/(kz)*Va*sum(krho(krho < 1).^2.*Bpar2_B02(krho<1))