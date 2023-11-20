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

clear;

list = {'0.01','0.02','0.05','0.1','0.2','0.3','0.4','0.6','0.8','1.0'};

%list = {'0.01','0.6','1.0'};
 
R_factor = [0.01;0.02;0.05;0.1;0.2;0.3;0.4;0.6;0.8;1.0];

strat_level = [19;14;13;14;14;13;12;18;12;18];

lengths_vs_month_up = zeros(10,1);

tic
for tt = 1:1

tt    
    
% Define Reduction Factor

r=R_factor(tt); 

% Define Constants

Omega = 7.2921e-5;
R = 6371000.0;
Ra = 287.04;
cp = 1005.7;
kappa = Ra/cp;
dt = 3600.0 * 6.0;
window_x = 1;
window_y = 1;

norm = 0;

path_temp = strcat('/net/aimsir/archive1/mkohl/interpolation/output_r_mode/data_rescale',list{tt},'/','temp_day0560h00.nc');
path_u = strcat('/net/aimsir/archive1/mkohl/interpolation/output_r_mode/data_rescale',list{tt},'/','ug_day0560h00.nc');
path_v = strcat('/net/aimsir/archive1/mkohl/interpolation/output_r_mode/data_rescale',list{tt},'/','vg_day0560h00.nc');
path_omega = strcat('/net/aimsir/archive1/mkohl/interpolation/output_r_mode/data_rescale',list{tt},'/','omega_day0560h00.nc');

% Read in Grid-Variables Lat/Lon 

lat_series_original = ncread(path_temp, 'latitude');
lon_series_original = ncread(path_temp, 'longitude');
latmax = 90; latmin = -90; lonmax = 360; lonmin = 0;
[lat_indices, lon_indices] = latlonindices(lat_series_original, lon_series_original, latmin, latmax, lonmin, lonmax);
plevels = double(ncread(path_temp, 'level'));


lon_series = lon_series_original(lon_indices);
lat_series = lat_series_original(lat_indices);
dphi = (lat_series(2) - lat_series(1)) / 180.0 * 3.1415926;
dlambda = (lon_series(2) - lon_series(1)) / 180.0 * 3.1415926;

% Specify Inversion Domain

event_latspan = [83:111]; %83 -111
event_lonspan = [1:256]; % 1:256
event_timespan = [191:194];
level_indices = [1:30]; % 5-30;
level = plevels(level_indices);
lat = lat_series_original(event_latspan);
lon = lon_series_original(event_lonspan);
phi = lat / 180.0 * 3.1415926;
lambda = lon / 180.0 * 3.1415926;
[~, Phi] = meshgrid(lambda, phi);

% Read in Variables (Lat/Lon need to be flipped)

start = [event_lonspan(1), event_latspan(1), level_indices(1), event_timespan(1)];
count = [length(event_lonspan), length(event_latspan), length(level_indices), length(event_timespan)];

omega_in= ncread(path_omega,'omega_plevel', start, count);

[T, u, v, omega,vor] = ...
        deal(zeros(length(event_latspan), length(event_lonspan), length(level_indices), length(event_timespan)));
 
for k = 1 : length(level_indices)
    for t = 1 : length(event_timespan)
        omega(:, :, k, t) = omega_in(:, :, k, t)';
    end
end



for t = 1:4



% find various levels

[~,p500] = min(abs(level-50000));
[~,p400] = min(abs(level-40000));
[~,p700] = min(abs(level-70000));
[~,plow] = min(abs(level-92500));
surf = 1;



norm = norm+1;

w = - omega(:,:,p500);

lengths = [];

for kk = 1:size(w,1)

w_slice = w(kk,:);
w_slice(w_slice>0) = 0;
w_down = w_slice<0;

input_array = w_down; 
input_array = [0 input_array 0];

input_diff = diff(input_array);

a0 = find(input_diff==1);
af = find(input_diff==-1);

dphi = abs(lon(2)-lon(1));
diffy = (af-a0)*2*pi*R*cosd(lat(kk))*dphi/360;

lengths = cat(2,lengths,diffy);
clear('input_array');

end

lengths_up = [];

for kk = 1 : size(w,1)

w_slice = w(kk,:);
w_slice(w_slice<0) = 0;
w_up = w_slice>0;

input_array = w_up; 
input_array = [0 input_array 0];

input_diff = diff(input_array);

a0 = find(input_diff==1);
af = find(input_diff==-1);

dphi = abs(lon(2)-lon(1));
diffy = (af-a0)*2*pi*R*cosd(lat(kk))*dphi/360;

lengths_up = cat(2,lengths_up,diffy);
clear('input_array');

end


end % end loop over time

toc


lengths_vs_month(tt) = mean(lengths);
lengths_vs_month_alternative(tt) = mean(1./lengths);
lengths_vs_month_up(tt) = mean(lengths_up);
lengths_vs_month_up_alternative(tt) = mean(1./lengths_up);

end

% convert index length to km
dphi = abs(lon(2)-lon(1));
%lengths = lengths*2*pi*R*cosd(45)*dphi/360;
%lengths_up = lengths_up*2*pi*R*cosd(45)*dphi/360;
lengths_vs_month = lengths_vs_month*2*pi*R*cosd(45)*dphi/360;
lengths_vs_month_up = lengths_vs_month_up*2*pi*R*cosd(45)*dphi/360;

lengths_vs_month_alternative = lengths_vs_month_alternative/(2*pi*R*cosd(45)*dphi/360);
lengths_vs_month_up_alternatve = lengths_vs_month_up_alternative/(2*pi*R*cosd(45)*dphi/360);


Ld = 5.8*1e5;

figure
histogram(Ld./lengths_up);