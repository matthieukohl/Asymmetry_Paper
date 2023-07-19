% Calculate the Asymmetry Based x average for ncep2;
% adapted from Paul's script

close all; clear;

years      = [2009:2018];
start_year = years(1);
end_year   = years(end);

months = {'jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov',...
    'dec'};

[lambda_ncep2_x_wzm] = deal(zeros(73,17,12));

for ii = 1:length(months)

season = months{ii};

norm = 0;
lambda = 0;

for year = years

disp('year:')
disp(year)

file = ['/net/aimsir/archive1/mkohl/era_interim/Asymmetry_Project/ncep2/omega.', num2str(year),'.nc'];


lat = ncread(file,'lat');
lon = ncread(file,'lon');
level = ncread(file,'level');
time = ncread(file,'time');

time_range = season_selector(file, season, start_year, end_year);

for time_counter = time_range

start = [1, 1, 1, time_counter];
count = [length(lon), length(lat), length(level), 1];

start_sf = [1, 1, time_counter];
count_sf = [length(lon), length(lat),1];

% have format lon/lat/level/time

omega = ncread(file,'omega',start,count);

% zonal-mean removed

omega = omega-mean(omega,1);

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

lambda_ncep2_x_wzm(:,:,ii) = lambda;


end

 save('data_ncep2/asymmetry_ncep2_x_wzm.mat', ...
      'lambda_ncep2_x_wzm', 'lat', 'level')

% [~,lat30] = min(abs(lat-40));
% [~,lat70] = min(abs(lat-60));
% 
% [~,p500] = min(abs(level-500));
% 
% lambda_paul = squeeze(mean(lambda_ncep2_x(lat70:lat30,p500,:),1));
% 
% figure(1)
% plot(lambda_paul,'Linewidth',1.5); hold on
% ylabel('\lambda')
% set(gca,'xtick',1:12,...
%  'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
% set(gca,'Fontsize',12)
% legend('Paul'); legend boxoff
% title('Lambda 30-70 NH')