function div_A = div(A_x, A_y, Phi, dphi, dlambda)

    R = 6371000.0;

    %div_A = 1 ./ (R * cos(Phi)) .* (d_dlambda(A_x, dlambda) + d_dphi(A_y .* cos(Phi), dphi));% + 
            %... %1 / R * tan(Phi) .* A_y + ...
            %1 / R * d_dphi(A_y, dphi);

% previous formula           
            
%     div_A = 1 ./ (R * cos(Phi)) .* d_dlambda(A_x, dlambda) + ...
%             1 / R * tan(Phi) .* A_y + ...
%             1 / R * d_dphi(A_y, dphi);

% should be minus sign in tan term 03/21

   div_A = 1 ./ (R * cos(Phi)) .* d_dlambda(A_x, dlambda)...
            + 1 ./ (R * cos(Phi)) .* d_dphi(A_y .* cos(Phi), dphi);

    return
