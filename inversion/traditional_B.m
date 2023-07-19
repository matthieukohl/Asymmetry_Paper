function B = traditional_B(phi, level, u, v, T, event_timespan, dphi, dlambda)

    B = zeros(size(u));
    Ra = 287.04;

    for t = 1 : length(event_timespan)
        for k = 1 : length(level)
            u_del_T = v_del(phi, u(:, :, k, t), v(:, :, k, t), T(:, :, k, t), dphi, dlambda);
            B(:, :, k, t) = Ra ./ level(k) .* spherical_laplacian(u_del_T, phi, dphi, dlambda);
        end
    end
end
   

