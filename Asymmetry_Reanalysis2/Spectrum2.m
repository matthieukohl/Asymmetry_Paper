function [Wavenumber_avg,spec] = Spectrum2(omega,lat,lon)

%omega(lat,lon,level)

W = zeros(size(omega,2),size(omega,3));
W_power_place = zeros(size(omega,2),size(omega,3));
Wavenumber = zeros(length(lat),1);
R = 6371000.0; % Planetary Radius 

for jj = 1:length(lat)

for k = 1:size(omega,3)
    
    W(:,k) = fft(omega(jj,:,k));
    W_power_place(:,k) = abs(W(:,k)).^2;
end

%W_power_place = squeeze(mean(W_power_place(:,:,:),2));


%% calculate Dominant Wavenumber 
W_power = mean(W_power_place(:,:),2);
Power(jj,:) = W_power;
R_eff = cosd(lat(jj))*R;
N = length(lon);

k = 2*pi/(2*pi*R_eff)*[0:N/2,-N/2+1:-1];

indx = N/2;

k_centroid = sum(k(1:indx)'.*W_power(1:indx))/sum(W_power(1:indx));

Wavenumber(jj) = k_centroid;

end

Wavenumber_avg = mean(Wavenumber);

spec = squeeze(mean(Power,1));

%Make Plot of Average Spectrum
% R_eff = cosd(50)*R; kk = 2*pi/(2*pi*R_eff)*(-N/2:N/2-1);
% 
% semilogx(R_eff*kk,squeeze(mean(Power,1))); hold on;

end