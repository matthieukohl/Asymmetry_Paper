function [T, ug, vg, omega, z, omega_b, tag] = event_read_Wrapper(Ny, Nx, Nz, ...
        event_year, input_path, event_latspan, event_lonspan, event_day, level_indices, string_1)

    [T, ug, vg, omega, z] = deal(nan(Ny, Nx, Nz, 3));
    omega_b = nan(Ny, Nx, Nz);
    if mod(event_day, 365) == 0.25
    
        time_indices = 365 / 0.25;
        input_year = event_year - 1;
        [T(:, :, :, 1), ug(:, :, :, 1), vg(:, :, :, 1), omega(:, :, :, 1), z(:, :, :, 1), omega_b] =...
                event_read_NetCDF_with_b(input_year, input_path, event_latspan, event_lonspan, time_indices, level_indices);
        time_indices = [1, 2];
        input_year = event_year;
        [T(:, :, :, 2:3), ug(:, :, :, 2:3), vg(:, :, :, 2:3), omega(:, :, :, 2:3), z(:, :, :, 2:3), ~] =...
                event_read_NetCDF_with_b(input_year, input_path, event_latspan, event_lonspan, time_indices, level_indices);
    
    elseif mod(event_day, 365) == 0
    
        time_indices = 365 / 0.25 + [-1, 0];
        input_year = event_year;
        [T(:, :, :, 1:2), ug(:, :, :, 1:2), vg(:, :, :, 1:2), omega(:, :, :, 1:2), z(:, :, :, 1:2), omega_b] =...
                event_read_NetCDF_with_b(input_year, input_path, event_latspan, event_lonspan, time_indices, level_indices);
        time_indices = 1;
        input_year = event_year + 1;
        if (input_year == 2000 && strcmp(string_1, 'historical')) || ...
           (input_year == 2100 && strcmp(string_1, 'rcp85'))
            disp(['year ', num2str(input_year), ' does not exist in the data, abandoned']);
            tag = false;
            return 
        end
        [T(:, :, :, 3), ug(:, :, :, 3), vg(:, :, :, 3), omega(:, :, :, 3), z(:, :, :, 3), ~] =...
        event_read_NetCDF_with_b(input_year, input_path, event_latspan, event_lonspan, time_indices, level_indices);
    
    else
    
        time_indices = mod(event_day, 365) / 0.25 + [-1, 0, 1];
        input_year = event_year;
        [T, ug, vg, omega, z, omega_b] = ...
                event_read_NetCDF_with_b(input_year, input_path, event_latspan, event_lonspan, time_indices, level_indices);
    
    end
    tag = true;

end
