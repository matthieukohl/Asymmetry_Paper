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

% res 1.5
if tf_hemisphere ==1
[lambda_era5_NH_x_coarse_grain] = deal(zeros(27,37,12));
else
[lambda_era5_SH_x_coarse_grain] = deal(zeros(27,37,12));
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

% % coarse grain omega: from 0.25 res to 1.25 res
% 
% omega = movmean(movmean(omega,6,1),6,2); % averaging
% 
% omega = omega(1:5:end,1:5:end,:); % downsampling

% coarse grain omega: from 0.25 res to 1.5 res

omega = movmean(movmean(omega,7,1),7,2); % averaging

omega = omega(1:6:end,1:6:end,:); % downsampling


% % coarse grain omega: from 0.25 res to 1.75 res
% 
% omega = movmean(movmean(omega,8,1),8,2); % averaging
% 
% omega = omega(1:7:end,1:7:end,:); % downsampling



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
lambda_era5_NH_x_coarse_grain(:,:,ii) = lambda;
else
lambda_era5_SH_x_coarse_grain(:,:,ii) = lambda; 
end

% lat downsampling
%lat = lat(1:10:end);
%lat = lat(1:5:end);
lat = lat(1:6:end);
%lat = lat(1:7:end);

end

if tf_hemisphere ==1

 save('data_era5/asymmetry_era5_NH_x_true_coarse_grain_1p5_res.mat', ...
      'lambda_era5_NH_x_coarse_grain', 'lat', 'level')
  
else
   save('data_era5/asymmetry_era5_SH_x_true_coarse_grain_1p5_res.mat', ...
      'lambda_era5_SH_x_coarse_grain', 'lat', 'level')  
end

% [~,lat30] = min(abs(lat-30));
% [~,lat70] = min(abs(lat-70));
% 
% [~,p500] = min(abs(level-500));
% 
% lambda_paul = squeeze(mean(lambda_erai_NH_x(lat70:lat30,p500,:),1));
% 
% figure(1)
% plot(lambda_paul,'Linewidth',1.5); hold on
% ylabel('\lambda')
% set(gca,'xtick',1:12,...
%  'xticklabel',{'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
% set(gca,'Fontsize',12)
% legend('Paul'); legend boxoff
% title('Lambda 30-70 NH')


% test fields
w_500 = omega(:,:,22);
w_500_coarse = movmean(movmean(w_500,11,1),11,2);
cmax = max(abs(w_500(:)));
nint = 20;
cint = cmax/nint;

figure
contourf(w_500',[-cmax:cint:cmax],'EdgeColor','none');
caxis([-cmax cmax]);
colormap(redblue(nint));
colorbar

figure
contourf(w_500_coarse',[-cmax:cint:cmax],'EdgeColor','none');
caxis([-cmax cmax]);
colormap(redblue(nint));
colorbar