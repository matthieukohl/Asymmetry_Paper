clear;
latmax = 90; % degrees
latmin = -90;
lonmax = 360; % degrees
lonmin = 0;
historicalpath = ['/net/chickpea/volume1/scratch/jgdwyer/cmip5/nomads.gfdl.noaa.gov/dods-data/CMIP5/output1/', ...
                        'NOAA-GFDL/GFDL-CM3/historical/3hr/atmos/3hr/r1i1p1/v20110601/pr/'];
historical_files = [historicalpath, 'pr_3hr_GFDL-CM3_historical_r1i1p1_1980010100-1984123123.nc'];
input_path = '/net/chickpea/volume1/scratch/ziweili/test1_GFDL/output_data/';

lat_series_original = ncread(historical_files, 'lat');
lon_series_original = ncread(historical_files, 'lon');
[lat_indices, lon_indices] = latlonindices(lat_series_original, lon_series_original, latmin, latmax, lonmin, lonmax);
plevels = double(ncread([input_path, 'z_historical_', num2str(1990), '.nc'], 'level'));

Omega = 7.2921e-5;
R = 6371000.0;
Ra = 287.04;
cp = 1005.7;
kappa = Ra/cp;
dt = 3600.0 * 6.0;
window_x = 1;
window_y = 1;

lon_series = lon_series_original(lon_indices);
lat_series = lat_series_original(lat_indices);
dphi = (lat_series(2) - lat_series(1)) / 180.0 * 3.1415926;
dlambda = (lon_series(2) - lon_series(1)) / 180.0 * 3.1415926;


input_path = '/net/chickpea/volume1/scratch/ziweili/test1_GFDL/output_data/';
input_year = 1990;
event_latspan = [20:36];
event_lonspan = [20:36];
event_timespan = [1:4];
%time_indices = [1:4];
level_indices = [5:30];
level = plevels(level_indices);
lat = lat_series_original(event_latspan);
lon = lon_series_original(event_lonspan);
phi = lat / 180.0 * 3.1415926;
lambda = lon / 180.0 * 3.1415926;
[~, Phi] = meshgrid(lambda, phi);
[T, ug, vg, omega, z, omega_b] =...
    event_read_NetCDF_with_b(input_year, input_path, event_latspan, event_lonspan, event_timespan, level_indices);

f0 = 2 * Omega * sin(lat_series_original(28) / 180 * 3.1415926);
beta = 2 * Omega / R * cos(lat_series_original(28) / 180 * 3.1415926);
[A, B, C, J, J1, J2, Qx, Qy, ...
    dug_dlambda, dvg_dlambda, dug_dphi, dvg_dphi, ...
    dT_dlambda, dT_dphi, sigma_accu] = deal(zeros(length(event_latspan), length(event_lonspan), length(level), length(event_timespan)));                  

[sigma, sigma_eff, dtheta_dp_ma, T_avg] = deal(zeros(length(level), length(event_timespan)));
[rhs, omega_QG] = deal(zeros(length(event_latspan), length(event_lonspan), length(level), length(event_timespan)));

for t = 1 : length(event_timespan)

    T_avg(:, t) = reshape(mean(mean(T(:, :, :, t), 1), 2), size(level));

    % effective static stability: sigma_eff
    temp_omega = squeeze(omega(:, :, :, t));
    [~, ~, theta_avg] = eff_stat_stab(level, T_avg(:, t), 1); % this routine is just used to calculate potential temperature                

                    % static stability: sigma
                    % the following formula is more accurate since it takes into account the moist part in the exponent 
                    % when calculating the potential temperature
                    sigma(:, t) =  - Ra * T_avg(:, t) ./ (level .* theta_avg) .* gradient(theta_avg, level);

                    % spatially varying static stability: sigma_accu
                    % this quantity is mainly used to calculate J term
                    Level = repmat(reshape(level, 1, 1, length(level)), length(event_latspan), length(event_lonspan), 1);
                    theta = T(:, :, :, t) .* (1e5 ./ Level) .^ kappa;
                    sigma_accu(:, :, :, t) = - Ra * T(:, :, :, t) ./ (Level .* theta) .* d_dp(theta, level);
                    
                    % deal with points where sigma_accu is not smooth
                    sigma_smooth_window = 3;    
                    if sigma_smooth_window ~= 0
                            for k = 1 : length(level)
                                sigma_accu_smoothed(:, :, k, t) = smooth2a(sigma_accu(:, :, k, t), sigma_smooth_window, sigma_smooth_window);
                            end
                        else
                            sigma_accu_smoothed(:, :, :, t) = sigma_accu(:, :, :, t);
                    end
                    % term A (Q vector)

                    for k = 1 : length(level)

                        % probably the ug and vg field need smoothing as well...
                        dug_dlambda(:, :, k, t) = d_dlambda(ug(:, :, k, t), dlambda);
                        dvg_dlambda(:, :, k, t) = d_dlambda(vg(:, :, k, t), dlambda);
                        dug_dphi(:, :, k, t) = d_dphi(ug(:, :, k, t), dphi);
                        dvg_dphi(:, :, k, t) = d_dphi(vg(:, :, k, t), dphi);

                        dT_dlambda(:, :, k, t) = d_dlambda(T(:, :, k, t), dlambda);
                        dT_dphi(:, :, k, t) = d_dphi(T(:, :, k, t), dphi);
                        Qx(:, :, k, t) = 1 / R^2 * (...
                            (1 ./ cos(Phi).^2) .* dug_dlambda(:, :, k, t) .* dT_dlambda(:, :, k, t) + ...
                            (1 ./ cos(Phi))    .* dvg_dlambda(:, :, k, t) .* dT_dphi(:, :, k, t));
                        Qy(:, :, k, t) = 1 / R^2 * (...
                            (1 ./ cos(Phi))    .* dug_dphi(:, :, k, t)    .* dT_dlambda(:, :, k, t) + ...
                                                  dvg_dphi(:, :, k, t)    .* dT_dphi(:, :, k, t));
                        div_temp = div(Qx(:, :, k, t), Qy(:, :, k, t), Phi, dphi, dlambda);
                        A(:, :, k, t) = 2 * Ra / level(k) * div_temp;
                    end
                    clear('temp_x', 'temp_y', 'div_temp');

                    % term B (planetary vorticity)

                    temp = d_dp(vg(:, :, :, t), level);
                    for j = 1 : length(lat)
                        B(j, :, :, t) = f0 * beta * temp(j, :, :);
                    end
                    clear('temp');

                    % term C (diabatic heating)

                    if t == 1
                        T1 = T(:, :, :, t);
                    end
                    if t < length(event_timespan)
                        T2 = T(:, :, :, t + 1);
                    end
                    for k = 1 : length(level)
                        % calculate diabatic heating using temperature field and geostrophic wind
                        %temp = v_del(phi, ug(:, :, k, t), vg(:, :, k, t), T(:, :, k, t), dphi, dlambda) - ...
                        %        sigma_accu(:, :, k, t) .* omega(:, :, k, t) / Ra * level(k);
                        %temp = temp * cp;
                                %sigma(k, t) * omega(:, :, k, t) / Ra * level(k);
                        if t == length(event_timespan)
                            temp = (T1(:, :, k) - T0(:, :, k)) / dt;
                        elseif t == 1
                            temp = (T2(:, :, k) - T1(:, :, k)) / dt;
                        else
                            temp = (T2(:, :, k) - T0(:, :, k)) / (2 * dt);
                        end

                        J1(:, :, k, t) = (temp + v_del(phi, ug(:, :, k, t), vg(:, :, k, t), T(:, :, k, t), dphi, dlambda)) * cp;
                        J2(:, :, k, t) = - sigma_accu(:, :, k, t) .* omega(:, :, k, t) / Ra * level(k) * cp;
                        J(:, :, k, t) = J1(:, :, k, t) + J2(:, :, k, t);
                        %J(:, :, k, t) = (J(:, :, k, t) + temp);

                        C(:, :, k, t) = - kappa / level(k) * spherical_laplacian(J(:, :, k, t), phi, dphi, dlambda);
                        %C1(:, :, k, t) = spherical_laplacian(sigma(k, t) * omega(:, :, k, t), phi, dphi, dlambda);
                        %C2(:, :, k, t) = C(:, :, k, t) - C1(:, :, k, t);

                    end

                    T0 = T1;
                    T1 = T2;

                    % right hand side

                    rhs(:, :, :, t) = A(:, :, :, t) + B(:, :, :, t) + C(:, :, :, t);


                                         % SIP_diagnostic(R, f0, Ni,                   Nj,            Nk, phi, lambda, level, ...
                                                           %sigma,       dphi, dlambda, rhs,  omega)
                    omega_QG(:, :, :, t) = SIP_diagnostic(R, f0, length(lon), length(phi), length(level), phi, lambda, level, ...
                                                            sigma(:, t), dphi, dlambda, rhs,  omega_b);
                    
                    %SIP_sigma(R, f0, length(lon), length(phi), length(level), phi, lambda, level, ...
                                                       % sigma_accu_smoothed(:, :, :, t), dphi, dlambda, rhs(:, :, :, t),  omega_b);
                end

figure
contour(omega_QG(:, :, 12, 2))
title('QG omega')

figure
contour(omega(:, :, 12, 2))
title('omega')

figure
contour(rhs(:, :, 12, 2))
title('right hand side')

