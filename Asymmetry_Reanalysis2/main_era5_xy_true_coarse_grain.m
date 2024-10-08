% Calculate the Asymmetry Based x average for erai;
% adapted from Paul's script

close all; clear;

years = [2009:2018];
%years = [2017];
start_year = years(1);
end_year   = years(end);

months = {'jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov',...
    'dec'};

hemisphere = 'south';
tf_hemisphere = strcmp(hemisphere,'north');

% res 1.25
if tf_hemisphere ==1
[lambda_era5_NH_xy_coarse_grain] = deal(zeros(37,12));
else
[lambda_era5_SH_xy_coarse_grain] = deal(zeros(37,12));
end

% % res 2.5
% if tf_hemisphere ==1
% [lambda_era5_NH_x_coarse_grain] = deal(zeros(17,37,12));
% else
% [lambda_era5_SH_x_coarse_grain] = deal(zeros(17,37,12));
% end

for ii = 1:length(months)

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

for time_counter = time_range

start = [1, 1, 1, time_counter];
count = [length(lon), length(lat), length(level), 1];

start_sf = [1, 1, time_counter];
count_sf = [length(lon), length(lat),1];

% have format lon/lat/level/time

omega = ncread(file,'w',start,count);

% % coarse grain omega: from 0.25 res to 2.5 res
% 
% omega = movmean(movmean(omega,11,1),11,2); % averaging
% 
% omega = omega(1:10:end,1:10:end,:); % downsampling

% coarse grain omega: from 0.25 res to 1.25 res

% omega = movmean(movmean(omega,6,1),6,2); % averaging
% 
% omega = omega(1:5:end,1:5:end,:); % downsampling
% 
% lat_sampled = lat(1:5:end);
% lon_sampled = lon(1:5:end);

% coarse grain omega: from 0.25 res to 1.25 res

omega = movmean(movmean(omega,7,1),7,2); % averaging

omega = omega(1:6:end,1:6:end,:); % downsampling

lat_sampled = lat(1:6:end);
lon_sampled = lon(1:6:end);

% % coarse grain omega: from 0.25 res to 1.75 res
% 
% omega = movmean(movmean(omega,8,1),8,2); % averaging
% 
% omega = omega(1:7:end,1:7:end,:); % downsampling
% 
% lat_sampled = lat(1:7:end);
% lon_sampled = lon(1:7:end);

omega_up = min(omega,0);

% Calculate Various Terms needed for Lambda

% co_var = mean(omega.*omega_up,1);
% omega_var = mean(omega.^2,1);
% omega_up_var = mean(omega_up.^2,1);
% 
% omega_avg = mean(omega,1);
% omega_up_avg = mean(omega_up,1);

% weighted

weight = repmat(cosd(lat_sampled),1,length(lon_sampled),length(level));
weight = reshape(weight,[length(lon_sampled),length(lat_sampled),length(level)]);

co_var = sum(sum(weight.*omega.*omega_up,1),2)./sum(sum(weight,1),2);
omega_var = sum(sum(weight.*omega.^2,1),2)./sum(sum(weight,1),2);
omega_up_var = sum(sum(weight.*omega_up.^2,1),2)./sum(sum(weight,1),2);

omega_avg = sum(sum(weight.*omega,1),2)./sum(sum(weight,1),2);
omega_up_avg = sum(sum(weight.*omega_up,1),2)./sum(sum(weight,1),2);

co_var_eddy = co_var-omega_up_avg.*omega_avg;
omega_var_eddy = omega_var-omega_avg.^2;

lambda = lambda + co_var_eddy./omega_var_eddy;

norm = norm+1;

%lambda_test = lambda_test+Lambda_x(omega);

end % end loop over time

end % end loop over years

lambda = squeeze(lambda/norm);

if tf_hemisphere ==1
lambda_era5_NH_xy_coarse_grain(:,ii) = lambda;
else
lambda_era5_SH_xy_coarse_grain(:,ii) = lambda; 
end

% lat downsampling
%lat = lat(1:10:end);
%lat = lat(1:5:end);
lat = lat(1:6:end);
%lat = lat(1:7:end);

end

if tf_hemisphere ==1

 save('data_era5/asymmetry_era5_NH_xy_true_coarse_grain_1p5_res.mat', ...
      'lambda_era5_NH_xy_coarse_grain', 'lat', 'level')
  
else
   save('data_era5/asymmetry_era5_SH_xy_true_coarse_grain_1p5_res.mat', ...
      'lambda_era5_SH_xy_coarse_grain', 'lat', 'level')  
end