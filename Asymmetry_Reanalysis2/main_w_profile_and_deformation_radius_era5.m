% Calculate the reduction factor and predicted asymmetry
% Based on x average for era5

close all; clear;

% Define Constants

Omega = 7.2921e-5; % Rotation Rate
R = 6371000.0; % Planetary Radius 
Ra = 287.04; % Gas-constant Air
Rc = 461.92; % Gas-constant Water Vapour
cp = 1005.7; % Heat Capacity Air
cpc = 1850; % Heat Capacity Water Vapour
kappa = Ra/cp; 
T_triple = 273.16; % Triple Point Temperature
p_triple = 611; % Triple Point Pressure
p0 = 101325; %US-Standard Atmosphere
L = 2493*10^3; % Latent Heat Water Vapour
g = 9.80665;  
epsilon = 0.621; % Mass ratio Water Vapour / Air

% Input for Era5 files

years  = [2009:2009];
%years = [2014:2014];
start_year = years(1);
end_year   = years(end);

months = {'jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov',...
    'dec'};

hemisphere = 'south';
tf_hemisphere = strcmp(hemisphere,'north');

if tf_hemisphere ==1
[w_up_NH,S_NH,Ld_NH] = deal(zeros(37,12));
else
[w_up_SH,S_SH,Ld_SH] = deal(zeros(37,12));
end

tic
for ii = 1:length(months)
ii
season = months{ii};

norm = 0;
counter_w_up = 0;
w_up = zeros(37,1);
S = zeros(37,1);

for year = years

disp('year:')
disp(year)

if tf_hemisphere ==1
file = ['/net/aimsir/archive1/mkohl/era_interim/Asymmetry_Project/era5/era5_NH_',num2str(ii),'_', num2str(year),'.nc'];
else
file = ['/net/aimsir/archive1/mkohl/era_interim/Asymmetry_Project/era5/era5_SH_',num2str(ii),'_',num2str(year),'.nc'];
end

lat = ncread(file,'latitude'); lat = double(lat); lat = flip(lat);
lon = ncread(file,'longitude');lon = double(lon);
level = ncread(file,'level'); level = double(level)*100; level = flip(level);

% Define f-factor
f = 2*Omega*sind(lat); % Coriolis parameter
beta = 2*Omega/R*cosd(lat); % beta-factor

time = ncread(file,'time');

time_range = season_selector(file, season, start_year, end_year);

tic

for t = time_range

start = [1, 1, 1, t];
count = [length(lon), length(lat), length(level), 1];

start_sf = [1, 1, t];
count_sf = [length(lon), length(lat),1];

% have format lon/lat/level/time

omega_in = ncread(file,'w',start,count);
T_in = ncread(file,'t',start,count);

T = rearrange_era(T_in,lat,lon,level);
omega = rearrange_era(omega_in,lat,lon,level);

clear('omega_in','T_in');

Level = repmat(reshape(level, 1, 1, length(level)), length(lat), length(lon), 1);
theta = T .* (p0 ./ Level).^ kappa;
sigma_accu = -Ra * T./ (Level .* theta) .* d_dp(theta, level);

r = r_parameter1(T,Level,level);
density = rho(T,Level);
lapse = -g*density.*d_dp(T,level)*10^3; % K/km

% find various levels

[~,p500] = min(abs(level-50000));
[~,p400] = min(abs(level-40000));
[~,p700] = min(abs(level-70000));
[~,plow] = min(abs(level-92500));
surf = 1;

% calculate the full vertical profile

S = S + squeeze(mean(mean(sigma_accu,1),2));

% calculate the vertical ascent profile

for tt = 1:size(omega,1)
 for jj = 1:size(omega,2)
     
 if any(omega(tt,jj,p500)<0)==1
     w_up = w_up + squeeze(omega(tt,jj,:));
     counter_w_up = counter_w_up + 1;
 end
 
  end
end


norm = norm+1;

end % end loop over time

toc

end % end loop over years

S = S/norm;
w_up = w_up/counter_w_up;


% calculate the deformation radius

w_up_pp = d_d2p(w_up,level);

w_pp_over_w = w_up_pp./w_up;

a = -mean(f.^2)./S;
a = a.*w_pp_over_w;

Ld = 1./sqrt(a);



if tf_hemisphere ==1
S_NH(:,ii) = S;    
w_up_NH(:,ii) = w_up;
Ld_NH(:,ii) = Ld;

else
S_SH(:,ii) = S;
w_up_SH(:,ii) = w_up;
Ld_SH(:,ii) = Ld;

end


end
toc

if tf_hemisphere ==1

 save('data_theory/deformation_radius_era5_NH_april_18th_500hPa.mat', ...
      'S_NH','w_up_NH','Ld_NH','lat','level')
  
else
   save('data_theory/deformation_radius_era5_SH_april_18th_500hPa.mat', ...
      'S_SH','w_up_SH','Ld_SH','lat','level')
    
end

% if tf_hemisphere ==1
% 
%  save('data_theory/deformation_radius_era5_NH_april_18th.mat', ...
%       'S_NH','w_up_NH','Ld_NH','lat','level')
%   
% else
%    save('data_theory/deformation_radius_era5_SH_april_18th.mat', ...
%       'S_SH','w_up_SH','Ld_SH','lat','level')
%     
% end

% make estimates of k_down

[a,p500] = min(abs(level-500*1e2));
w = - omega(:,:,p500);
points_wup = w(w<0); % number of indices with ascent at 500hPA
total_points = size(w,1)*size(w,2); % total number of indices in the domain

area_up_fraction = size(points_wup,1)/total_points; % estimate area fraction
total_area = areaquad(lat(1),lon(1),lat(end),lon(end)); % calculate quadrangle of unit sphere
total_area = total_area * 4*pi*(R/1000)^2; % in km^2

area_up = area_up_fraction*total_area;


% measure length of index

lengths = [];

for ii = 1: size(w,1)

w_slice = w(ii,:);
w_slice(w_slice>0) = 0;
w_down = w_slice<0;

input_array = w_down; 
input_array = [0 input_array 0];

input_diff = diff(input_array);

a0 = find(input_diff==1);
af = find(input_diff==-1);

lengths = cat(2,lengths,af-a0);
clear('input_array');

end

% convert index length to km
dphi = abs(lon(2)-lon(1));
lengths = lengths*2*pi*R/1000*cosd(45)*dphi;

figure
histogram(lengths)

Ld = 1e6/(2*sqrt(2));
lengths = lengths/Ld;

figure
histogram(lengths)


% identify length of strings

% Define input array
input_array = [1 0 1 0 1 0 1 0 1 1 0 1 1 1 0 0 1 1 1 1];
input_array = [0 input_array 0];

input_diff = diff(input_array);

a0 = find(input_diff==1);
af = find(input_diff==-1);

lengths = af-a0;

histogram(lengths)


        
     
