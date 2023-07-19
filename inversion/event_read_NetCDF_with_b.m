function [T_out, ug_out, vg_out, omega_out, z_out, omega_b_out] = ...
    event_read_NetCDF_with_b(input_year, input_path, event_latspan, event_lonspan, time_indices, level_indices)

    g = 9.8;

    if input_year < 2000
        string_1 = 'historical';
    elseif input_year >= 2080
        string_1 = 'rcp85';
    else
        disp(['year ', num2str(input_year), ' does not exist in the data, abandoned']);
        return;   
    end
    
    start = [1, event_latspan(1), level_indices(1), time_indices(1)];
    count = [Inf, length(event_latspan), length(level_indices), length(time_indices)];

    T = ncread([input_path, 'ta_', string_1, '_', num2str(input_year), '.nc'], 'ta_plevel', start, count);
    T = T(event_lonspan, :, :, :);
    ug = ncread([input_path, 'ug_', string_1, '_', num2str(input_year), '.nc'], 'ug', start, count);
    ug = ug(event_lonspan, :, :, :);
    vg = ncread([input_path, 'vg_', string_1, '_', num2str(input_year), '.nc'], 'vg', start, count);
    vg = vg(event_lonspan, :, :, :);
    omega = ncread([input_path, 'omega_', string_1, '_', num2str(input_year), '.nc'], 'omega_plevel', start, count);
    omega = omega(event_lonspan, :, :, :);
    z = ncread([input_path, 'z_', string_1, '_', num2str(input_year), '.nc'], 'z_plevel', start, count);
    z = z(event_lonspan, :, :, :) * g; % this is important because the geopotential calculated does not have factor g

    % add climatologically averaged omega as lateral boundary condition
    omega_path = '/net/chickpea/volume1/scratch/ziweili/test1_GFDL/omega/';
    start = [1, event_latspan(1), level_indices(1)];
    count = [Inf, length(event_latspan), length(level_indices)];
    omega_b = ncread([omega_path, 'omega_avg_', string_1, '.nc'], 'omega_avg', start, count);
    omega_b = omega_b(event_lonspan, :, :);
    
    [T_out, ug_out, vg_out, omega_out, z_out] = ...
        deal(zeros(length(event_latspan), length(event_lonspan), length(level_indices), length(time_indices)));
    omega_z_out = zeros(length(event_latspan), length(event_lonspan), length(level_indices));
    for k = 1 : length(level_indices)
        for t = 1 : length(time_indices)
            T_out(:, :, k, t) = T(:, :, k, t)';
            ug_out(:, :, k, t) = ug(:, :, k, t)';
            vg_out(:, :, k, t) = vg(:, :, k, t)';
            omega_out(:, :, k, t) = omega(:, :, k, t)';
            z_out(:, :, k, t) = z(:, :, k, t)';
        end
        omega_b_out(:, :, k) = omega_b(:, :, k)';
    end


end
