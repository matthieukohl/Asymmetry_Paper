function [Wavenumber] = Spectrum(omega,lat,lon)

%W = zeros(size(omega,2),size(omega,3),size(omega,4));
%W_power_place = zeros(size(omega,2),size(omega,3),size(omega,4));

W = zeros(size(omega,2),1,size(omega,4));
W_power_place = zeros(size(omega,2),1,size(omega,4));
R = 6371000.0; % Planetary Radius 

indx = find(lat == 46);

for t = 1 : size(omega,4)
  
for k = 1:1
    
    W(:,k,t) = fftshift(fft(omega(indx,:,16,t)));
    W_power_place(:,k,t) = abs(W(:,k,t)).^2;
end

W_power_place = squeeze(mean(W_power_place(:,:,:),2));

end


%% calculate Dominant Wavenumber 
W_power = mean(W_power_place(:,:),2);

R_eff = cosd(50)*R;
N = length(lon);
k = 2*pi/(2*pi*R_eff) *(-N/2:N/2-1);

k_centroid = sum(k(90:end).*W_power(90:end)')/sum(W_power(90:end));

Wavenumber = R_eff*k_centroid;

end

