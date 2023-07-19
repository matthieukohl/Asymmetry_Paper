% Calculate the Asymmetry Based x average for era5 - coarse
% resolution

close all; clear;

years      = [2009:2018];
start_year = years(1);
end_year   = years(end);

months = {'jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov',...
    'dec'};

hemisphere = 'south';
tf_hemisphere = strcmp(hemisphere,'north');

if tf_hemisphere ==1
[lambda_era5_NH_x_coarse] = deal(zeros(17,12));
else
[lambda_era5_SH_x_coarse] = deal(zeros(17,12));
end

for ii = 1:length(months)

season = months{ii};

norm = 0;
lambda = 0;

for year = years

disp('year:')
disp(year)

if tf_hemisphere ==1
file = ['/net/aimsir/archive1/mkohl/era_interim/Asymmetry_Project/era5/omega_2009_2018_NH_coarse_resolution.nc'];
else
file = ['/net/aimsir/archive1/mkohl/era_interim/Asymmetry_Project/era5/omega_2009_2018_SH_coarse_resolution.nc'];
end

lat = ncread(file,'latitude');
lon = ncread(file,'longitude');
time = ncread(file,'time');

time_range = season_selector(file, season, start_year, end_year);

for time_counter = time_range

start = [1, 1, time_counter];
count = [length(lon), length(lat),1];

start_sf = [1, 1, time_counter];
count_sf = [length(lon), length(lat),1];

% have format lon/lat/level/time

omega = ncread(file,'w',start,count);

omega_up = min(omega,0);

% Calculate Various Terms needed for Lambda

co_var = mean(omega.*omega_up,1);
omega_var = mean(omega.^2,1);
omega_up_var = mean(omega_up.^2,1);

omega_avg = mean(omega,1);
omega_up_avg = mean(omega_up,1);

co_var_eddy = co_var-omega_up_avg.*omega_avg;
omega_var_eddy = omega_var-omega_avg.^2;

lambda = lambda + co_var_eddy./omega_var_eddy;

norm = norm+1;

%lambda_test = lambda_test+Lambda_x(omega);

end % end loop over time

end % end loop over years

lambda = squeeze(lambda/norm);

if tf_hemisphere ==1
lambda_era5_NH_x_coarse(:,ii) = lambda;
else
lambda_era5_SH_x_coarse(:,ii) = lambda; 
end


end

if tf_hemisphere ==1

 save('data_era5/asymmetry_era5_NH_x_coarse.mat', ...
      'lambda_era5_NH_x_coarse','lat')
  
else
   save('data_era5/asymmetry_era5_SH_x_coarse.mat', ...
      'lambda_era5_SH_x_coarse','lat')  
end