
% make sample plot of vertical velocity field

close all; clear;

%years      = [2009:2018];
years = [2017];
start_year = years(1);
end_year   = years(end);

months = {'jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov',...
    'dec'};

hemisphere = 'north';
tf_hemisphere = strcmp(hemisphere,'north');

if tf_hemisphere ==1
[lambda_era5_NH_x] = deal(zeros(161,37,12));
else
[lambda_era5_SH_x] = deal(zeros(161,37,12));
end

for ii = 6:6

season = months{ii};

norm = 0;
lambda = 0;

for year = years

disp('year:')
disp(year)

if tf_hemisphere ==1
file = ['/net/aimsir/archive1/mkohl/era_interim/Asymmetry_Project/era5/era5_NH_',num2str(ii),'_', num2str(year),'.nc'];
else
file = ['/net/aimsir/archive1/mkohl/era_interim/Asymmetry_Project/era5/era5_SH_',num2str(ii),'_',num2str(year),'.nc'];
end

lat = ncread(file,'latitude');
lon = ncread(file,'longitude');
level = ncread(file,'level');
time = ncread(file,'time');

time_range = season_selector(file, season, start_year, end_year);

for time_counter = 2

start = [1, 1, 1, time_counter];
count = [length(lon), length(lat), length(level), 1];

start_sf = [1, 1, time_counter];
count_sf = [length(lon), length(lat),1];

% have format lon/lat/level/time

omega = ncread(file,'w',start,count);

end
end
end

% lat lon plot

[Lon,Lat] = meshgrid(lon,lat);

figure('Renderer', 'painters', 'Position', [10 10 1000 300]); 
w = -omega(:,:,22);
cmax = max(abs(w(:))); 
nint = 20; cint = cmax/nint; 
contourf(Lon,Lat,w',-cmax:cint:cmax,'EdgeColor','none'); 
colormap(redblue(nint)); colorbar; 
caxis([-cmax cmax]);
xlabel('Lon({\circ})')
ylabel('Lat({\circ})')
title('\rm Vertical Velocity (500hPa)');
set(gca,'FontSize',12);

% make a 2D power spectrum plot

% Load the SHTOOLS library
addpath('shtools')

% Define the degree and order of the spherical harmonics
lmax = 10;

% Define the grid of longitude and latitude points
nlon = 360;
nlat = 180;
lon = linspace(0, 360, nlon+1)';
lat = linspace(-90, 90, nlat+1)';
[lon, lat] = meshgrid(lon(1:end-1), lat(1:end-1));
lon = lon*pi/180;
lat = lat*pi/180;

% Define the latitude range of interest
lat_start = lat(1);
lat_end = lat(end);

% Find the indices of the latitude range of interest
lat_idx = find(lat >= lat_start*pi/180 & lat <= lat_end*pi/180);

% Load the variable f(lon,lat,pressure)
f = omega;

% Select the latitudes of interest from the variable
f_lat = f(:,lat_idx,:);

% Calculate the power spectrum for each pressure level
power_spectrum = zeros(lmax+1,1);
for i = 1:size(f_lat,3)
    % Apply the spherical harmonic transform to the data
    [clm, ~, ~] = shtridap(squeeze(f_lat(:,:,i))', lmax, 'irr', 'n', length(lat_idx));
    % Calculate the power spectrum
    power_spectrum_i = shpower(clm);
    % Add the power spectrum to the total power spectrum
    power_spectrum = power_spectrum + power_spectrum_i;
end

% Sum the power spectrum over all pressure levels to obtain the total power spectrum
total_power_spectrum = sum(power_spectrum, 2);

% Plot the power spectrum
l = 0:lmax;
figure
loglog(l, total_power_spectrum)
xlabel('Degree l')
ylabel('Power spectrum')





% slice plot

indxl = 100;
indxr = 700;
[a,indy] = min(abs(lat-45));
[a,indp] = min(abs(level-500));

figure('Renderer', 'painters', 'Position', [10 10 500 300])

plot(lon(indxl:indxr),-omega(indxl:indxr,indy,indp),'b','Linewidth',1.4); hold on; 
plot(lon(indxl:indxl),zeros(indxr-indxl,1),'k:'); 
xlim([lon(indxl) lon(indxr)]);
ylim([-1.5 1.5])
xlabel('Longitude (degrees)'); 
ylabel('-\omega (Pa/s)'); box off; 
set(gca,'FontSize',12)
saveas(gcf,'/disk7/mkohl/Mode_Kohl_OGorman/figures_sls/w_slice_era5','epsc')
