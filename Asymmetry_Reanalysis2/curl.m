function rot_A = curl(phi, Ax, Ay, dphi, dlambda)

    r = 6371000.0;
    sizeA = size(Ax);
    Ni = sizeA(2);
    Nj = sizeA(1);
    Omega = 7.2921e-5;
    [dAy_dx, dAx_dy] = deal(zeros(sizeA));

    % x (longitude) boundary condition
    dAy_dx(:, 1)  = 1 ./ (r * cos(phi(:))) .* ( - 3/2 * Ay(:, 1)  + 2 * Ay(:, 2)      - 1/2 * Ay(:, 3))      / dlambda;
    dAy_dx(:, Ni) = 1 ./ (r * cos(phi(:))) .* (   3/2 * Ay(:, Ni) - 2 * Ay(:, Ni - 1) + 1/2 * Ay(:, Ni - 2)) / dlambda;

    % y (latitude) boundary condition
    dAx_dy(1, :)  = 1 ./ (r * cos(phi(1))) .* ( - 3/2 * Ax(1, :) *cos(phi(1)) + 2 * Ax(2, :)     *cos(phi(2))      - 1/2 * Ax(3, :)     *cos(phi(3)))      / dphi;
    dAx_dy(Nj, :) = 1 ./ (r * cos(phi(Nj))).* (   3/2 * Ax(Nj, :)*cos(phi(Nj))- 2 * Ax(Nj - 1, :)*cos(phi(Nj - 1)) + 1/2 * Ax(Nj - 2, :)*cos(phi(Nj - 2))) / dphi; 
    
    % differentiate in x direction, vg and delf_x
    
    for i = 2 : Ni - 1
        for j = 1 : Nj
            dAy_dx(j, i) = dAy_dx(j, i) + 1 / (r * cos(phi(j))) * (Ay(j, i + 1) - Ay(j, i - 1)) / (2 * dlambda);
        end
    end
    
    
    for i = 1 : Ni
        for j = 2 : Nj - 1
            dAx_dy(j, i) = dAx_dy(j, i) + 1 / (r * cos(phi(j))) * (Ax(j + 1, i)*cos(phi(j + 1)) - Ax(j - 1, i)*cos(phi(j - 1))) / (2 * dphi);
        end
    end
    
    rot_A = dAy_dx - dAx_dy;

    return
    

