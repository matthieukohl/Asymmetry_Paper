function dAx_dy = d_dy(phi, Ax, dphi)

    r = 6371000.0;
    sizeA = size(Ax);
    Ni = sizeA(2);
    Nj = sizeA(1);
    [dAx_dy] = deal(zeros(sizeA));

    % y (latitude) boundary condition
    dAx_dy(1, :)  = 1 ./ (r * cos(phi(1))) .* ( - 3/2 * Ax(1, :) *cos(phi(1)) + 2 * Ax(2, :)     *cos(phi(2))      - 1/2 * Ax(3, :)     *cos(phi(3)))      / dphi;
    dAx_dy(Nj, :) = 1 ./ (r * cos(phi(Nj))).* (   3/2 * Ax(Nj, :)*cos(phi(Nj))- 2 * Ax(Nj - 1, :)*cos(phi(Nj - 1)) + 1/2 * Ax(Nj - 2, :)*cos(phi(Nj - 2))) / dphi; 
    
    % differentiate in y-direction
   
    for i = 1 : Ni
        for j = 2 : Nj - 1
            dAx_dy(j, i) = dAx_dy(j, i) + 1 / (r * cos(phi(j))) * (Ax(j + 1, i)*cos(phi(j + 1)) - Ax(j - 1, i)*cos(phi(j - 1))) / (2 * dphi);
        end
    end
    

return