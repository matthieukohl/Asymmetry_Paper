function dAy_dx = d_dx(phi, Ay, dlambda)

    r = 6371000.0;
    sizeA = size(Ay);
    Ni = sizeA(2);
    Nj = sizeA(1);
    Omega = 7.2921e-5;
    [dAy_dx] = deal(zeros(sizeA));

    % x (longitude) boundary condition
    dAy_dx(:, 1)  = 1 ./ (r * cos(phi(:))) .* ( - 3/2 * Ay(:, 1)  + 2 * Ay(:, 2)      - 1/2 * Ay(:, 3))      / dlambda;
    dAy_dx(:, Ni) = 1 ./ (r * cos(phi(:))) .* (   3/2 * Ay(:, Ni) - 2 * Ay(:, Ni - 1) + 1/2 * Ay(:, Ni - 2)) / dlambda;
    
    % differentiate in x direction
    
    for i = 2 : Ni - 1
        for j = 1 : Nj
            dAy_dx(j, i) = dAy_dx(j, i) + 1 / (r * cos(phi(j))) * (Ay(j, i + 1) - Ay(j, i - 1)) / (2 * dlambda);
        end
    end
    
return