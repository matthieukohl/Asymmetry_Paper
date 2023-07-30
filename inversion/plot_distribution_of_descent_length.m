% Inverts the 3D-Omega equation del^2(r(w)N^2w)+f^2w_zz = div Q + beta * dv/dz
% for a given RHS (r-simulations) iteratively using Dirichlet-BC w=0 
% on the boundaries.

clear;

list = {'_true0.01','0.02','0.05','0.1','0.2','0.3','0.4','_true0.6','0.8','_true1.0'};
list_end = {'0.01','0.02','0.05','0.1','0.2','0.3','0.4','0.6','0.8','1.0'};

%list = {'0.01','0.6','1.0'};
 
r_gcm = [0.01;0.02;0.05;0.1;0.2;0.3;0.4;0.6;0.8;1.0];

strat_level = [19;14;13;14;14;13;12;18;12;18];

w_up_vs_r = zeros(30,length(r_gcm));
w_up_vs_r_normalized = zeros(30,length(r_gcm));
Ld_vs_r = zeros(30,length(r_gcm));

tic
for tt = 1:1

tt        

% Define Constants

Omega = 7.2921e-5;
R = 6371000.0;
Ra = 287.04;
cp = 1005.7;
kappa = Ra/cp;
dt = 3600.0 * 6.0;
window_x = 1;
window_y = 1;

% Deformation radius

Ld_new = zeros(30,1);

sigma_new = zeros(30,1);

w_up = zeros(30,1);

counter = 0;
counter_w_up = 0;

for i= 9:9 %9:12
i
% Define Path 

path_temp = strcat('/disk7/mkohl/interpolation/output1/data_rescale',list{tt},'/temp_day0',num2str(25 * i, '%.3d'),'h00.nc');
path_u = strcat('/disk7/mkohl/interpolation/output1/data_rescale',list{tt},'/ug_day0',num2str(25 * i, '%.3d'),'h00.nc');
path_v = strcat('/disk7/mkohl/interpolation/output1/data_rescale',list{tt},'/vg_day0',num2str(25 * i, '%.3d'),'h00.nc');
path_omega = strcat('/disk7/mkohl/interpolation/output1/data_rescale',list{tt},'/omega_day0',num2str(25 * i, '%.3d'),'h00.nc');

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
event_timespan = [1:100];
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

T_in= ncread(path_temp,'temp_plevel', start, count);
u_in= ncread(path_u,'ug', start, count);
v_in= ncread(path_v,'vg', start, count);
omega_in= ncread(path_omega,'omega_plevel', start, count);

[T, u, v, omega,omega_pp,vor] = ...
        deal(zeros(length(event_latspan), length(event_lonspan), length(level_indices), length(event_timespan)));
 
for k = 1 : length(level_indices)
    for t = 1 : length(event_timespan)
        T(:, :, k, t) = T_in(:, :, k, t)';
        u(:, :, k, t) = u_in(:, :, k, t)';
        v(:, :, k, t) = v_in(:, :, k, t)';
        omega(:, :, k, t) = omega_in(:, :, k, t)';
    end
end

% start calculating

[a,p500] = min(abs(level-500e2));
    
w = - squeeze(omega(:,:,p500,:));

lengths = [];

for t = 1:size(omega,4)

for ii = 1: size(w,1)

w_slice = squeeze(w(ii,:,t));
w_slice(w_slice>0) = 0;
w_down = w_slice<0;

input_array = w_down; 
input_array = [0 input_array 0];

input_diff = diff(input_array);

a0 = find(input_diff==1);
af = find(input_diff==-1);

lengths = cat(2,lengths,af-a0);
clear('input_array','a0','af');

end

end

lengths_up = [];

for t = 1:size(omega,4)

for ii = 1: size(w,1)

w_slice = w(ii,:,t);
w_slice(w_slice<0) = 0;
w_up = w_slice>0;

input_array = w_up; 
input_array = [0 input_array 0];

input_diff = diff(input_array);

a0 = find(input_diff==1);
af = find(input_diff==-1);

lengths_up = cat(2,lengths_up,af-a0);
clear('input_array');

end

end



end

end

distance_index = 2*pi*R*cosd(45)*(lon(2)-lon(1))/360;
lengths = lengths*distance_index;
lengths_up = lengths_up*distance_index;
%Ld = 1e6/(2*sqrt(2));

Ld = 5.79*1e5; % using S500 from model; see file plot_w_spec_r_simulation

lengths = lengths/Ld;
lengths_up = lengths_up/Ld;

figure('Renderer', 'painters', 'Position', [10 10 1200 400])

subplot(1,2,1)

histogram(lengths)
set(gca,'Yscale','log');
xlabel('Ld')
set(gca,'FontSize',12)
title(['Ld=',num2str(round(mean(lengths),2)),' kd=',num2str(round(1/mean(lengths),2)),' (',num2str(round(mean(1./lengths),2)),') r=',num2str(r_gcm(tt))])

subplot(1,2,2)

histogram(lengths_up)
set(gca,'Yscale','log');
xlabel('Lu')
set(gca,'FontSize',12)
title(['Lu=',num2str(round(mean(lengths_up),2)),' ku=',num2str(round(1/mean(lengths_up),2)),' (',num2str(round(mean(1./lengths_up),2)),') r=',num2str(r_gcm(tt))])

ku = 1/mean(lengths_up);
kd = 1/mean(lengths);

ku = mean(1./lengths_up); 
kd = mean(1./lengths);

k = 2*ku*kd/(ku+kd);


lengths = 1./lengths;

figure
histogram(lengths)
set(gca,'Yscale','log');
xlabel('Ld')
set(gca,'FontSize',12)
title(['kd=',num2str(mean(lengths))])


counter = counter + 1;








