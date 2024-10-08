% Calculate the Asymmetry Based x average for ncep2;
% adapted from Paul's script

close all; clear;

years      = [2009:2018];
start_year = years(1);
end_year   = years(end);

months = {'jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov',...
    'dec'};

[lambda_ncep2_xy] = deal(zeros(17,12));

hemisphere = 'south';
tf_hemisphere = strcmp(hemisphere,'north');

for ii = 1:length(months)

season = months{ii};

norm = 0;
lambda = 0;

for year = years

disp('year:')
disp(year)

file = ['/net/aimsir/archive1/mkohl/era_interim/Asymmetry_Project/ncep2/omega.', num2str(year),'.nc'];


lat = ncread(file,'lat');

if tf_hemisphere==1
    
   [a,lat30] = min(abs(lat-30));
   [a,lat70] = min(abs(lat-70));
   lat = lat(lat70:lat30);
else
   [a,lat30] = min(abs(lat+30));
   [a,lat70] = min(abs(lat+70));
   lat = lat(lat30:lat70);
end
    

lon = ncread(file,'lon');
level = ncread(file,'level');
time = ncread(file,'time');

time_range = season_selector(file, season, start_year, end_year);

for time_counter = time_range
    
    if tf_hemisphere==1

    start = [1, lat70, 1, time_counter];
    count = [length(lon), lat30-lat70+1, length(level), 1];

    else
        
        
    start = [1, lat30, 1, time_counter];
    count = [length(lon), lat70-lat30+1, length(level), 1];
  
    end

start_sf = [1, 1, time_counter];
count_sf = [length(lon), length(lat),1];

% have format lon/lat/level/time

omega = ncread(file,'omega',start,count);

omega_up = min(omega,0);

% Calculate Various Terms needed for Lambda

% weighted

weight = repmat(cosd(lat),1,length(lon),length(level));
weight = reshape(weight,[length(lon),length(lat),length(level)]);

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

lambda_ncep2_xy(:,ii) = lambda;


end

if tf_hemisphere ==1
save('data_ncep2/asymmetry_ncep2_xy_NH.mat', ...
      'lambda_ncep2_xy', 'lat', 'level')
  
else
    
save('data_ncep2/asymmetry_ncep2_xy_SH.mat', ...
      'lambda_ncep2_xy', 'lat', 'level')
end
    