function D = Deformation(phi, level, u, v, T, event_timespan, dphi, dlambda)
    
[ux,uy,vx,vy,Tx,Ty] = deal(zeros(size(u)));


D = zeros(size(u));
Ra = 287.04;

for t = 1 : length(event_timespan)
    for k = 1 : length(level)
    ux(:,:,k,t) = d_dx(phi, u(:,:,k,t), dlambda);
    uy(:,:,k,t) = d_dy(phi, u(:,:,k,t), dphi);
    vx(:,:,k,t) = d_dx(phi, v(:,:,k,t), dlambda);
    vy(:,:,k,t) = d_dy(phi, v(:,:,k,t), dphi);
    Tx(:,:,k,t) = d_dx(phi, T(:,:,k,t), dlambda);
    Ty(:,:,k,t) = d_dy(phi, T(:,:,k,t), dphi);
    end
end

for t = 1 : length(event_timespan)
    for k = 1 : length(level)
D(:,:,k,t) = 2*Ra ./ level(k).*(v_del(phi, ux(:, :, k, t), uy(:, :, k, t), Tx(:, :, k, t), dphi, dlambda)+...
    v_del(phi, vx(:, :, k, t), vy(:, :, k, t), Ty(:, :, k, t), dphi, dlambda));      
    end
end
    
    
end
   